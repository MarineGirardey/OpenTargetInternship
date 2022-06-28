import requests
import pyspark.sql.functions as f
from pyspark.sql import SparkSession
import argparse
from json import JSONDecodeError
import logging
import re
import pandas as pd
import sys


def main():

    spark = SparkSession.builder.getOrCreate()

    # MOLECULES
    molecule_df = (
        spark.read
        .parquet(args.parquet)
        .select(
            f.col('inchiKey').alias('inchikey'), f.col('id').alias('chemblId'), 
            'name', 'linkedTargets', 'linkedDiseases'
            )
    )
    logging.info(f'STATISTICS -> Number in molecules : {molecule_df.count()}')

    # INCHIKEY MOLECULES
    inchikey_df = (
        spark.read
            .csv(args.inchikey, sep=',', header=True, comment='#')
            .select(
                f.col('InChIKey').alias('inchikey'), 
                f.col('CCD_ID').alias('pdbCompoundId')
                )
    )
    logging.info(f'STATISTICS -> Number in inchikeys : {inchikey_df.count()}')

    # MOLECULE WITH COMPOUND ID
    molecules_inchikey_join = (
        molecule_df
        .join(inchikey_df, on='inchikey')
        # .sample(0.001, 3)
    )
    logging.info(f'STATISTICS -> Number of our molecules in inchikeys : {molecules_inchikey_join.count()}')

    # MOLECULE WITH PDB STRUCTURES
    molecules_w_pdb = (
        spark
        .createDataFrame(molecules_inchikey_join
                         .toPandas()
                         .assign(pdbStructureId=lambda df: df.pdbCompoundId.apply(get_structure))
                         )
        .filter(f.size('pdbStructureId') != 0)
    )
    logging.info(f'STATISTICS -> Number of molecules with one or more pdb structure : {molecules_w_pdb.count()}')

    # NON PERTINENT COMPOUNDS
    excluded_compounds = (
        spark
        .createDataFrame(
            pd.read_csv('https://zhanggroup.org/BioLiP/ligand_list', header=None)
            .rename(columns=({0: 'pdbCompoundId'}))
            .astype(str)
        )
        .persist()
    )

    # NON PERTINENT COMPOUNDS EXCLUDED
    molecules_w_pdb_drug_filtered = (
        molecules_w_pdb
        .join(excluded_compounds, on='pdbCompoundId', how='left_anti')
    )
    logging.info(f'STATISTICS -> Number after removing non relevant compounds : {molecules_w_pdb_drug_filtered.count()}')

    # TARGET
    pdb_target_df = (
        spark.read.csv(args.ensembl, sep=',', header=True, comment='#')
        .select(
            f.col('PDB').alias('pdbStructureId'), 
            f.col('GENE_ID').alias('geneId'),
            f.col('CHAIN').alias('chainId'),
            f.col('SP_PRIMARY').alias('uniprot')
        )
        .distinct()
        .persist()
    )
    logging.info(f'STATISTICS -> Number of (ensembl) targets : {pdb_target_df.count()}')

    # COMPOUND STRUCTURE AND TARGET
    gene_mapped_structures = (
        molecules_w_pdb_drug_filtered
        .select(
            'pdbCompoundId',
            f.explode(f.col('pdbStructureId')).alias('pdbStructureId')
        )
        .groupby('pdbStructureId')
        .agg(f.collect_set(f.col('pdbCompoundId')).alias('compoundIds'))

        .join(pdb_target_df, on='pdbStructureId', how='left')
        .filter(f.col('geneId').rlike('ENSG\d+'))

        # Reorganizing the dataset:
        .groupby('pdbStructureId')
        .agg(
            f.collect_set(
                f.struct(
                    f.col('geneId'),
                    f.col('chainId'),
                    f.col('uniprot')
                )
            ).alias('chains'),
            f.first(f.col('compoundIds')).alias('compoundIds')
        )
        .persist()
    )
    logging.info(f'STATISTICS -> Number of structure with compounds and targets: {gene_mapped_structures.count()}')


    # CREATE JSON OUTPUT FILE
    (
        gene_mapped_structures
        .write.mode('overwrite').json(args.output + 'structure_drug_target_o.json')
    )

    if args.folder:
        molecules_inchikey_join.write.mode('overwrite').json(args.output + 'intermediary_data/molecules_inchikey_join.json')
        molecules_w_pdb.write.mode('overwrite').json(args.output + 'intermediary_data/molecules_w_pdb.json')
        molecules_w_pdb_drug_filtered.write.mode('overwrite').json(args.output + 'intermediary_data/molecules_w_pdb_drug_filtered.json')


def get_structure(pdb_compound_id: str) -> list:
    """Fetching structure identifiers from PDBkb REST API
    Args:
        pdb_compound_id: string, a single compound identifier
    Returns:
        List of PDB structure identifiers where the compound can be found
    """

    # URL of the API to obtain the structure id list
    url = f'https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/{pdb_compound_id}'
    response = requests.get(url)

    # Try to obtain the API return in a json format
    try:
        data = response.json()
        return data[pdb_compound_id]

    # Impossible to obtain the json format because there is an error
    except JSONDecodeError:
        print(f'Failed to return structures to: {pdb_compound_id}')
        if len(response.json()) == 0:
            return []

    # Key error the strucutre id is not good or not referenced
    except KeyError:
        logging.warning(f'NO DATA -> Empty data was returned for: {pdb_compound_id}')
        return []


if __name__ == '__main__':

    global spark

    program_description = '''
    This script retrieves the structures corresponding to the complex between molecules from parquet files and targets.
    How to find the files?
    molecule file here (Drug): http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.11/output/etl/parquet/molecule/*
    inchikey file : https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/components_inchikeys.csv
    ensembl file : https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_ensembl.csv
    '''

    parser = argparse.ArgumentParser(add_help=True, description=program_description)

    parser.add_argument('-p',
                        '--parquet',
                        help='Path to parquet folder',
                        default=None,
                        metavar='parquet_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-i',
                        '--inchikey',
                        help='Path to inchikey file',
                        default=None,
                        metavar='inchikey_file_path',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output',
                        help='Path to output folder',
                        default=None,
                        metavar='output_folder_path',
                        type=str,
                        required=True)

    parser.add_argument('-e',
                        '--ensembl',
                        help='Path to pdb_chain_ensembl.csv file',
                        default=None,
                        metavar='ensembl_file_path',
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
                        '--folder',
                        help='folder with intermediate data',
                        action='store_true',
                        required=False)

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=args.log_file,
        force=True)

    logging.StreamHandler(sys.stderr)

    logging.info(program_description)

    main()