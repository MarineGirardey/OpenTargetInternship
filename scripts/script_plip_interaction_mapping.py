import json
from json import JSONDecodeError

import requests
from functools import reduce
import pandas as pd

from plip.basic import config

from plip.structure.preparation import PDBComplex, PLInteraction
from plip.exchange.report import BindingSiteReport
from plip.basic import config

import dask.dataframe as dd
from pandarallel import pandarallel

from itertools import chain

import logging

import argparse

import sys
from dask.distributed import Client
import signal
import time


def main():

    logging.warning('Program begin.')

    client = Client()  # start distributed scheduler locally.  Launch dashboard

    config.DNARECEPTOR = True

    logging.warning('Loading input.')

    # Dataset witht all the details, produced earlier:
    input_dataset = (
        pd.read_csv(args.input_file, sep=",")
        .groupby("pdbStructureId")
        .agg(pd.unique)
    )

    # PARALLELISATION METHODS

    # ddf = dd.from_pandas(input_dataset, npartitions=args.nb_partitions)

    # input_dataset = (
    #     ddf
    #    .assign(
    #         new_col = ddf.map_partitions(
    #            lambda df: df.apply(lambda row: characerize_complex(row), axis=1), meta=(None, 'f8')
    #        )
    #        .map_partitions(lambda df: df.apply(run_plip, axis=1), meta=(None, 'f8'))
    #    )
    #    .compute(scheduler='processes')
    # )
    
    pandarallel.initialize()

    input_dataset['new_col'] = input_dataset.parallel_apply(
        characerize_complex, axis=1
    )

    logging.info('PLIP interaction computations finished.')

    plip_output_pd_df = pd.DataFrame(list(chain.from_iterable(
        input_dataset
        .loc[lambda df: df.new_col.apply(lambda x: len(x) >0)]
        .assign(new_col = lambda df: df.new_col.apply(lambda l: [value for value in l if value != {}]))
        .new_col
        .to_list()
    )))

    plip_output_pd_df.to_csv(args.output_file, index=False, header=True)

    logging.info('Program finished, file saved as "interaction_structure_drug_plip_output.csv".')


class GetPDB:
    
    PDB_URL = 'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent'
    
    def __init__(self, data_folder: str) -> None:
        self.data_folder = data_folder
        
    
    def get_pdb(self, pdb_structure_id: str) -> str:
        """Reading file from a given loaction fetch and save if not found"""
    
        logging.warning(f'Start searching the file.')
    
        try:
            # Readind data from the given location:
            with open(f'{self.data_folder}/pdb{pdb_structure_id}.ent', 'rt') as f:
                logging.warning(f'try.')
                data = f.read()
                logging.warning(f'try ok.')
    
        except FileNotFoundError:
            logging.warning(f'except.')
            # Fetch data from the web
            data = self.fetch_pdb(pdb_structure_id)
            
            # Save file
            with open(f'{self.data_folder}/pdb{pdb_structure_id}.ent', 'wt') as f:
                f.write(data)
                logging.warning(f'except ok.')

        logging.warning(f'File obtained.')
    
        return data
    

    def fetch_pdb(self, pdb_structure_id: str)-> str:
        """This function fetches the pdb file from ePDB server as a string

        Args:
            pdb_structure_id (str)
        Returns:
            structure data in pdb format as string eg 'AIN:A:1202'
        """
        data = ""
    
        headers={'Content-Type': 'text/plain'}
    
        if not pdb_structure_id:
            return ''

        try:
            response = requests.get(self.PDB_URL.format(pdb_structure_id), headers=headers)
            if response.headers['Content-Type'] != 'text/plain; charset=UTF-8':
                pass
            else:
                data = response.text
    
        except:
            data = ''
    
        logging.warning(f'File obtained or exception while scrapping.')

        return data


def parse_interaction(interaction: PLInteraction, compound_id:str, pdb_id:str) -> dict:

    interaction_type = interaction.__doc__.split('(')[0]
    
    if interaction_type == 'waterbridge':
        return {}

    # Parsing data form the interaction:
    return {
        'pdb_structure_id': pdb_id,
        'compound_id': compound_id,
        'interaction_type': interaction_type,
        'prot_residue_number': interaction.resnr,
        'prot_residue_type': interaction.restype,
        'prot_chain_id': interaction.reschain
    }



def characerize_complex(row):

    compounds = row[0]
    pdb_id = row.name

    # Get pdb data:
    # pdb_id = row['pdbStructureId']
    # compounds = row['pdbCompoundId']

    logging.warning(f'Characerize_complex: {pdb_id, compounds}')
    
    gpdb = GetPDB(data_folder=args.pdb_folder)

    pdb_data = gpdb.get_pdb(pdb_id)

    if pdb_data:

        # Load into plip:
        mol_complex = PDBComplex()
        
        try:
            mol_complex.load_pdb(pdb_data, as_string=True)

        except:
            pass
        
        if mol_complex.ligands:
            
            # Filtering out only the relevant ligands:
            ligands_of_interest = [ligand for ligand in mol_complex.ligands if ligand.hetid in compounds]
                
            # Characterizing relevant complex:
            [mol_complex.characterize_complex(ligand) for ligand in ligands_of_interest]
            
            logging.warning(f'Sites computation finished.')

            # Extract details from ligands:
            return [parse_interaction(interaction, compound.split(':')[0], pdb_id) for compound, interaction_set in mol_complex.interaction_sets.items() for interaction in interaction_set.all_itypes]

        else:
            return []
    
    else:
        return []


if __name__ == '__main__':

    program_description = '''
    Compute PLIP interactions between PDB structures and compounds (Drugs).
    '''

    parser = argparse.ArgumentParser(add_help=True, description=program_description)

    parser.add_argument('-i',
                        '--input_file',
                        help='Path to the csv file with structures and drugs to compute PLIP interactions on.',
                        default=None,
                        metavar='csv_structure_drug_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-o',
                        '--output_file',
                        help='Path to the output csv file with PLIP interactions computed for each structure-drug combinations.',
                        default=None,
                        metavar='csv_plip_output_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-f',
                        '--pdb_folder',
                        help='Path to the pdb folder with pdb files downloaded from pdbe website and used by PLIP to compute interactions.',
                        default=None,
                        metavar='pdb_folder_path',
                        type=str,
                        required=True)
    
    parser.add_argument('-p',
                        '--nb_partitions',
                        help='Number of Dask partitions (I build 30 partitions with 8 cores).',
                        default=None,
                        metavar='nb_dask_partitions',
                        type=int,
                        required=False)

    args = parser.parse_args()

    logging.warning(program_description)

    main()
