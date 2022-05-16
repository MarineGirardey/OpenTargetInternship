import requests
import argparse
import pandas as pd
import time
from functools import reduce
from pandarallel import pandarallel



def get_pdb_sifts_mapping(pdb_id: str) -> pd.DataFrame:
    URL = f'https://www.ebi.ac.uk/pdbe/graph-api/mappings/ensembl/{pdb_id}'
    
    while True:
        try:
            data = requests.get(URL).json()
            break
        except:
            time.sleep(1)
            continue
            
    # Let's try to process data. And shape it into a dataframe:
    try:
        return (
            pd.DataFrame(reduce(lambda x,y: x + y['mappings'], data[pdb_id]['Ensembl'].values(), []))
            .assign(
                author_start = lambda df: df.start.apply(lambda start: start['author_residue_number']),
                author_end = lambda df: df.end.apply(lambda end: end['author_residue_number']),
                uniprot_position = lambda df: df.apply(lambda row: list(range(row['unp_start'], row['unp_end']+1)), axis=1),
                diff = lambda df: df.apply(lambda row: row['author_start'] - row['unp_start'], axis=1)
            )
            .explode('uniprot_position')
            .assign(
                prot_residue_number = lambda df: df.apply(lambda row: row['uniprot_position'] + row['diff'], axis=1)
            )
            [['accession', 'chain_id', 'uniprot_position', 'prot_residue_number']]
            .rename(columns={'chain_id': 'prot_chain_id'})
            .drop_duplicates()
        )

    # If anything goes wrong we are returning and empty dataframe with the right columns:
    except:
        return pd.DataFrame(columns=['accession', 'uniprot_position', 'prot_chain_id', 'prot_residue_number'])


def map2uniprot(plip_df: pd.DataFrame) -> pd.DataFrame:
    # Extracting pdb identifier:
    pdb_id = plip_df.pdb_structure_id.iloc[0]
    
    # Fetch mappings from pdb api:
    sifts_df = get_pdb_sifts_mapping(pdb_id)
    
    # Join with mapping:
    return (
        plip_df
        .merge(sifts_df, on=['prot_chain_id', 'prot_residue_number'], how='left')
    )


def main(plip_file: str, id_mapping_file: str, output_file: str) -> None:

    # Reading plip output:
    df = pd.read_json(plip_file, orient='records', lines=True)
    print(f'Number of rows in the plip dataset: {len(df)}')
    print(f'Number of pdb structures in the plip dataset: {len(df.pdb_structure_id.unique())}')

    # Grouping data by pdb structure id:
    grouped = df.groupby('pdb_structure_id')
    # grouped = (
    #     df
    #     .query('pdb_structure_id == "3e7g" or pdb_structure_id == "13gs" ')
    #     .groupby('pdb_structure_id')
    # )
    print(f'Number of groups in the grouped dataset: {len(grouped)}')

    # Selecting one of the groups:
    pandarallel.initialize(progress_bar=True)
    mapped_df = grouped.parallel_apply(map2uniprot).reset_index(drop=True)
    print(f'Number of structures in the mapped dataset: {mapped_df.pdb_structure_id.unique()}')

    # Joining with translation id:
    uniprot2ensembl = (
        pd.read_csv(id_mapping_file, sep='\t', names=['accession', 'source', 'identifier'])
        .query('source == "Ensembl_PRO"')
        .drop('source', axis=1)
        .rename(columns={'identifier': 'translation_id'})
    )

    mapped_w_ensp = (
        mapped_df
        .merge(uniprot2ensembl, on='accession', how='inner')
    )

    mapped_w_ensp.to_json(output_file, compression='infer', orient='records', lines=True)
    print('Done.')



if __name__ == '__main__':
    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(
        description='This script maps plip output to Ensembl translations.')
    parser.add_argument('--output_file', '-o', type=str,
                        help='Output file. gzipped JSON', required=True)
    parser.add_argument('--plip_file', '-i', type=str,
                        help='A JSON file describing exeriment metadata', required=True)
    parser.add_argument('--id_mapping_file', type=str,
                        help='File with uniprot to ensembl translation id mappings', required=False)
    args = parser.parse_args()

    main(args.plip_file, args.id_mapping_file, args.output_file)