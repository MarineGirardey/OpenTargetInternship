import requests
import argparse
import pandas as pd
import time
from functools import reduce
from pandarallel import pandarallel
import os
import logging

import pyspark.sql.functions as f
from pyspark.sql import SparkSession
from pyspark.sql.types import StructField, StructType, IntegerType, StringType

spark = (
    SparkSession.builder
    .master('local[*]')
    .config("spark.driver.memory", "15g")
    .getOrCreate()
)


def get_pdb_sifts_mapping(pdb_id: str) -> pd.DataFrame:
    """Obtain conventional (uniprot) residue number from PDBe GRAPH API
    Args:
        pdb_id: string, pdb structure identifier
    Returns:
        DataFrame for each group of the groupby with conventional numbers of all residue in a structure
    """

    URL = f'https://www.ebi.ac.uk/pdbe/api/mappings/ensembl/{pdb_id}'

    # To obtain the mapping file even in case of bad connexion
    # TODO: Need to be improved because infinite loop
    while True:
        try:
            data = requests.get(URL).json()
            break
        # TODO: Bad exception
        except:
            time.sleep(1)
            continue

    # Process the data for each structure (grouby) mapping file. And shape it into a dataframe
    try:
        return (
            pd.DataFrame(reduce(lambda x, y: x + y['mappings'], data[pdb_id]['Ensembl'].values(), []))

            .assign(
                author_start=lambda df: df.start.apply(lambda start: start['author_residue_number']),
                author_end=lambda df: df.end.apply(lambda end: end['author_residue_number']),
                uniprot_position=lambda df: df.apply(lambda row: list(range(
                    row['unp_start'], row['unp_end'] + 1
                )), axis=1),
                diff=lambda df: df.apply(lambda row: row['author_start'] - row['unp_start'], axis=1)
            )

            .explode('uniprot_position')
            .assign(
                prot_residue_number=lambda df: df.apply(lambda row: row['uniprot_position'] + row['diff'], axis=1)
            )
            [['accession', 'chain_id', 'uniprot_position', 'prot_residue_number']]
            .rename(columns={'chain_id': 'prot_chain_id'})
            .drop_duplicates()
        )

    # If anything goes wrong we are returning and empty dataframe with the right columns:
    # TODO: Bad exception
    except:
        return pd.DataFrame(columns=['accession', 'uniprot_position', 'prot_chain_id', 'prot_residue_number'])


def map2uniprot(plip_df: pd.DataFrame) -> pd.DataFrame:
    """Intermediar function to lauch the 'get_pdb_sifts_mapping' for each structure (groupby)
    function and return one single dataframe
    Args:
        plip_df: pandas dataframe, dataframe in a groupby format to execute a function on each group
    Returns:
        one final DataFrame with all residues numbers with the conventional number of each structures
    """

    # Extracting pdb identifier:
    pdb_id = plip_df.pdb_structure_id.iloc[0]

    # Fetch mappings from pdb api:
    sifts_df = get_pdb_sifts_mapping(pdb_id)

    # Join with mapping:
    return (
        plip_df
        .merge(sifts_df, on=['prot_chain_id', 'prot_residue_number'], how='left')
    )


def main(plip_file: str, prot_mapping_file: str, output_file: str, loc_mapping_file: str) -> None:

    # Reading PLIP output
    df = (
        spark.read.json(plip_file)
        .toPandas()
    )

    logging.info(f'Number of rows in the plip dataset: {len(df)}')
    logging.info(f'Number of pdb structures in the plip dataset: {len(df.pdb_structure_id.unique())}')

    # Grouping data by pdb structure id:
    grouped = df.groupby('pdb_structure_id')

    logging.info(f'Number of groups in the grouped dataset: {len(grouped)}')

    # Selecting one of the groups:
    pandarallel.initialize(progress_bar=True)
    mapped_df = grouped.parallel_apply(map2uniprot).reset_index(drop=True)

    logging.info(f'Number of structures in the mapped dataset: {mapped_df.pdb_structure_id.unique()}')

    # Open the translation id mapping file (give the good translation id because there is error in the graph api)
    uniprot2ensembl = (
        pd.read_csv(prot_mapping_file, sep='\t', names=['accession', 'source', 'identifier'])
        .query('source == "Ensembl_PRO"')
        .drop('source', axis=1)
        .rename(columns={'identifier': 'translation_id'})
    )

    # Merge translation id df to replace the wrong translation id by the correct ones
    mapped_w_ensp = (
        mapped_df
        .merge(uniprot2ensembl, on='accession', how='inner')
    )

    # Just convert the pandas df to a pyspark df and put good column name
    mapped_w_ensp_sp = (
        spark.createDataFrame(mapped_w_ensp)
        .withColumnRenamed("pdb_structure_id", "pdbStructureId")
        .withColumnRenamed("prot_chain_id", "chainId")
        .withColumnRenamed("compound_id", "pdbCompoundId")
    )

    # Open mapping file for residue and genomic location that daniel generated with the genome anotated from gencode
    generated_mapping = (
        spark.read.csv(args.loc_mapping_file, sep="\t", header=True)
        .withColumnRenamed("protein_id", "translation_id")
        .withColumnRenamed("gene_id", "geneId")
        .withColumnRenamed("amino_acid_position", "prot_residue_number")
    )

    # Join df with conventional (uniprot) residue numbers and Daniel's mapping file
    residue_gen_pos = (
        mapped_w_ensp_sp
        .join(generated_mapping, on=['translation_id', 'prot_residue_number'], how='inner')
        .drop("uniprot_position")
        .toPandas()
    )

    # Save the final df with genomic location of the interacting residues
    residue_gen_pos.to_json(output_file, compression='infer', orient='records', lines=True)
    logging.info('Done.')


if __name__ == '__main__':

    program_description = '''
    Compute genomic location of interacting residues identified with PLIP in the previous script.
    '''

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(
        description='This script maps plip output to Ensembl translations.')

    parser.add_argument('--output_file', '-o', type=str,
                        help='Output file. gzipped JSON', required=True)

    parser.add_argument('--plip_file', '-i', type=str,
                        help='A csv file with plip output', required=True)

    parser.add_argument('--prot_mapping_file', '-p', type=str,
                        help='File with uniprot to ensembl translation id mappings', required=False)

    parser.add_argument('--loc_mapping_file', '-m', type=str,
                        help='Daniel mapping file for residue location and genomic location', required=True)

    parser.add_argument('-l',
                        '--log_file',
                        help='File to save the logs to.',
                        default=None,
                        metavar='log_file',
                        type=str,
                        required=True)

    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, force=True)
    logging.info(program_description)

    main(args.plip_file, args.prot_mapping_file, args.output_file, args.loc_mapping_file)
