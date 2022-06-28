import requests
import pandas as pd
import psutil
import logging
import argparse
import time

from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction
from plip.exchange.report import BindingSiteReport

from pyspark.sql import SparkSession
import pyspark.sql.functions as f

from pandarallel import pandarallel
from itertools import chain

import json
from json import JSONDecodeError

def main():

    # Important PLIP parameter I forgot why but fail if False
    config.DNARECEPTOR = True

    spark = (
        SparkSession.builder
        .master('local[*]')
        .config("spark.driver.memory", "15g")
        .getOrCreate()
    )

    logging.info('Loading input.')

    # Open input which is the output of the previous script (get structure + target + drug info)
    input_dataset = (
        spark.read.json(args.input_file)
        .select("pdbStructureId", "chains", "compoundIds")
        .toPandas()
    )

    # Pandarallel initialisation required
    pandarallel.initialize(
        nb_workers=psutil.cpu_count(),
        progress_bar=True,
    )

    # Use .apply from pandas combined with pandarallel because pyspark is not working with PLIP 
    # (object PLIP incompatible with pyspark)
    # 'new_col' will contain a list with interacting residue info from plip computations
    input_dataset['new_col'] = input_dataset.parallel_apply(
        characerize_complex, axis=1
    )

    logging.info('PLIP interaction computations finished.')

    # Kirill block code to put in a nice form the 'new_col', ask him for details
    plip_output_pd_df = pd.DataFrame(list(chain.from_iterable(
        input_dataset
        .loc[lambda df: df.new_col.apply(lambda x: len(x) > 0)]
        .assign(new_col=lambda df: df.new_col.apply(lambda l: [value for value in l if value != {}]))
        .new_col
        .to_list()
    )))

    # Groupby and organise to avoid repetition (many type of interaction for a same residue-residue contact)
    plip_output_df = (
        spark.createDataFrame(plip_output_pd_df)
        .groupby('pdb_structure_id', 'compound_id', 'prot_residue_number', 'prot_residue_type', 'prot_chain_id')
        .agg(
            f.collect_set(
                f.col('interaction_type')
            )
            .alias('interaction_type')
        )
    )

    plip_output_df.write.json(args.output_file)

    logging.info('!!!Program finished!!!')


class GetPDB:
    """Class to get the PDB file, download locally if not already in the pdb folder
    Objects:
        data_folder: path of the folder which contain the pdb structure files
    Returns:
        data: the PDB structure file opened and returned as a string
    """

    # URL to download the pdb files
    PDB_URL = 'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent'

    def __init__(self, data_folder: str) -> None:
        self.data_folder = data_folder

    def get_pdb(self, pdb_structure_id: str) -> str:
        """
        Reading file (save in a string) from a given loaction fetch and save if not found from API
        Objects:
            pdb_structure_id: id of the structure to get the PDB file
        Returns:
            data: the PDB structure file opened and returned as a string
        """
        logging.info('Start searching the file.')

        # TODO: We can try to NOT download the file in local but to get the string (might be quicker)
        try:
            # Readind data from the given LOCAL location
            with open(f'{self.data_folder}/pdb{pdb_structure_id}.ent', 'rt') as f:
                data = f.read()

        # If the file is not already downloaded
        except FileNotFoundError:

            # Fetch data from the web
            data = self.fetch_pdb(pdb_structure_id)
            # Save file
            with open(f'{self.data_folder}/pdb{pdb_structure_id}.ent', 'wt') as f:
                f.write(data)

        logging.info('File obtained.')

        return data

    def fetch_pdb(self, pdb_structure_id: str) -> str:
        """This function fetches the pdb file from ePDB server as a string
        Args:
            pdb_structure_id: PDB id of the structure
        Returns:
            data: pdb structure file as string
        """

        data = ""

        # Set format of the API answer to avoir error response
        headers = {'Content-Type': 'text/plain'}

        if not pdb_structure_id:
            return ''

        # Try to obtain the structure file from the API and save it in a variable as a string
        try:
            response = requests.get(self.PDB_URL.format(pdb_structure_id), headers=headers)
            if response.headers['Content-Type'] != 'text/plain; charset=UTF-8':
                pass
            else:
                data = response.text

        # TODO: Bad exception not enough specific
        # Here we want to avoid any error from the API response. commexion error 
        # or other error which return a text/plain text
        except:
            data = ''

        return data


def parse_interaction(interaction: PLInteraction, compound_id: str, pdb_id: str) -> dict:
    """Function to put in a nice format PLIP results
    Args:
        interaction: PLIP interaction object
        compound_id: drug id
        pdb_id: pdb structure id
    Returns:
        a dictionary to format result column with the interacting residue for each structure
    """

    # Get interaction tyoe from the interaction PLIP object
    interaction_type = interaction.__doc__.split('(')[0]

    # TODO: Ask Daniel why this filter
    if interaction_type == 'waterbridge':
        return {}

    # Parsing data form the interaction
    return {
        'pdb_structure_id': pdb_id,
        'compound_id': compound_id,
        'interaction_type': interaction_type,
        'prot_residue_number': interaction.resnr,
        'prot_residue_type': interaction.restype,
        'prot_chain_id': interaction.reschain
    }


def characerize_complex(row: tuple):
    """Function to run other functions to obtain the interacting residue from PLIP
    Args:
        row: line of the input dataframe which is the drug, structure, target info
    Returns:
        a dataframe with the interacting residue of each drug-target complex of the dataframe
    """
    # Start the timer
    start_time = time.time()

    pdb_id = row[0]
    compounds = row[2]
    result = []

    # Initialise and run the class to get pdb files of the corresponding structure
    gpdb = GetPDB(data_folder=args.pdb_folder)
    pdb_data = gpdb.get_pdb(pdb_id)

    # If we managed to obtain the file (not all structure can be obtained)
    if pdb_data:

        # Load into PLIP, initialise PLIP objets
        mol_complex = PDBComplex()

        try:
            # Load the pdb file as a string in PLIP objet
            mol_complex.load_pdb(pdb_data, as_string=True)

        # TODO: Avoid this, need to be change
        except:
            pass

        # We had cases where the complex/ligands are not found or something
        # TODO: Why we loose data here?
        if mol_complex.ligands:

            # From all ligands found in the structure, only keep the ligand of interest here
            ligands_of_interest = [ligand for ligand in mol_complex.ligands if ligand.hetid in compounds]

            # Characterizing drug-ligand complex (compute interactions)
            [mol_complex.characterize_complex(ligand) for ligand in ligands_of_interest]

            logging.info('Interaction computations finished.')

            # Extract details from ligands, technical but right
            result = [
                parse_interaction(
                    interaction, compound.split(':')[0], pdb_id
                )
                for compound, interaction_set in mol_complex.interaction_sets.items()
                for interaction in interaction_set.all_itypes
            ]

    # Log the complex info data: PDB ID, list of ligands, size of data, processing time in seconds.
    logging.info(f'COMPLEX_INFO\t{pdb_id}\t{compounds}\t{len(pdb_data)}\t{time.time()-start_time}')

    return result


if __name__ == '__main__':

    global spark

    program_description = '''
    Compute PLIP interactions between target and compounds from PDB structures.
    '''

    parser = argparse.ArgumentParser(add_help=True, description=program_description)

    parser.add_argument('-i',
                        '--input_file',
                        help='Path to the json file output of the first script (get structure, target infos from drugs)',
                        default=None,
                        metavar='json_structure_drug_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-o',
                        '--output_file',
                        help='Path to the output compressed json file with PLIP interactions computed for each structure-drug combinations.',
                        default=None,
                        metavar='json_plip_output_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-l',
                        '--log_file',
                        help='File to save the logs to.',
                        default=None,
                        metavar='log_file',
                        type=str,
                        required=True)

    parser.add_argument('-f',
                        '--pdb_folder',
                        help='Path to the pdb folder with pdb files downloaded from pdbe website and used by PLIP to compute interactions.',
                        default=None,
                        metavar='pdb_folder_path',
                        type=str,
                        required=True)

    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, force=True)

    logging.info(program_description)

    main()
