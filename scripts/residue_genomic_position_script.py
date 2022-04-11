import pyspark.sql.functions as f
from pyspark.sql import SparkSession
import pandas as pd
import findspark
findspark.init()
import pyspark # Call this only after findspark.init()
from pyspark.context import SparkContext

import requests
from json import JSONDecodeError
import logging
import argparse
from pandarallel import pandarallel
import psutil

sc = SparkContext.getOrCreate()
spark = SparkSession(sc)


def main():

    # # Progress bar removed because : OverflowError: int too big to convert
    # pandarallel.initialize(
    #     nb_workers=psutil.cpu_count(),
    #     progress_bar=True,
    # )

    # PLIP INPUT (contain wanted data)
    plip_json_input = (

        # "gene_mapped_structures.json"
        spark.read.json(args.plip_input)

        .select("pdbStructureId", "chains", f.explode("compoundIds").alias("pdbCompoundId"))

        .select("pdbStructureId", "pdbCompoundId", f.explode("chains").alias("chains"))

        )

    plip_json_input_v2 = (plip_json_input

                        .withColumn("chainId", plip_json_input["chains.chainId"])

                        .withColumn("geneId", plip_json_input["chains.geneId"])

                        .withColumn("uniprotId", plip_json_input["chains.uniprot"])

                        .drop("chains")

                        )

    # PLIP OUTPUT
    plip_csv_output = (

        # "output.csv"
        spark.read.csv(args.plip_output, header=True, sep=",")

        .withColumnRenamed("pdb_structure_id", "pdbStructureId")

        .withColumnRenamed("compound_id", "pdbCompoundId")

        .withColumnRenamed("prot_chain_id", "chainId")

        )

    # JOIN to have gene id (used for filter mapping file on the gene id)
    plip_output_target_id = (

        plip_json_input_v2

        .join(plip_csv_output, on=["pdbStructureId", "chainId", "pdbCompoundId"])

        .withColumnRenamed("interaction_type", "intType")

        .withColumnRenamed("prot_residue_number", "protResNb")

        .withColumnRenamed("prot_residue_type", "protResType")

    )

    # Aggregate
    plip_output_agg = (

        plip_output_target_id

        .groupby([f.col('geneId'),
                f.col('uniprotId'),
                f.col("pdbStructureId").alias("pdbStructId")
                ])

        .agg(f.collect_set(f.col("pdbCompoundId")).alias("pdbCompId"),

            f.collect_set(f.struct(
                f.col('intType'),
                f.col('chainId'),
                f.col('protResType'),
                f.col('protResNb')))

            .alias("intType, chain, resType, resNb")
            )
        )

    # Test set
    if args.test_set:

        plip_output_agg = (
            plip_output_agg
            .select("geneId", "pdbStructId", "intType, chain, resType, resNb")
            .filter(plip_output_agg.pdbStructId.rlike('1dqa'))

        )

    pandarallel.initialize(
        nb_workers=psutil.cpu_count(),
        progress_bar=True,
    )

    # Pandas Apply
    plip_output_agg_pd = plip_output_agg.toPandas()
    plip_output_agg_pd['res_infos'] = plip_output_agg_pd.parallel_apply(
        fetch_gapi_ensembl_mapping, axis=1
    )

    # Final DF
    final_df = (
        spark.createDataFrame(plip_output_agg_pd)
        .drop("intType, chain, resType, resNb")
        .select("geneId", "pdbStructId", f.explode("res_infos").alias("res_infos"))
    )

    final_df.toPandas().to_csv(args.output_folder + "/residue_genomic_position.csv", index=False, header=True)


def fetch_gapi_ensembl_mapping(row):
    """This function fetches the graph api ensembl mapping file from ePDB server

    Args:
        rows of the DataFrame (each structures)
    Returns:
        a column for each structure with genomic positions and other infos about residues
    """

    gene_id = row[0]
    pdb_struct_id = row[2]
    residue_info = row[4]

    url = f'https://www.ebi.ac.uk/pdbe/graph-api/mappings/ensembl/{pdb_struct_id}'
    headers={'Content-Type': 'application/json'}

    try:
        response = requests.get(url, headers=headers)

        if response.headers['Content-Type'] != 'application/json':
            e_mapping_file = None

        else:
            e_mapping_file = response.json()

    except JSONDecodeError:
        if len(response.json()) == 0:
            e_mapping_file = None

    except KeyError:
        e_mapping_file = None

    if e_mapping_file:
        return filter_dict_file(e_mapping_file, pdb_struct_id, gene_id, residue_info)
    else:
        return None


def filter_dict_file(e_mapping_file, pdb_struct_id, gene_id, residue_info):
    """This function compute genomic positions of each residues of a structure

    Args:
        e_mapping_file : to extract genomic position of a range a residues
        pdb_struct_id : the concerned structure
        gene_id : gene ensembl id
        residue_info : all residues to extract genomic position and info about them
    Returns:
        all_residu__info_list : List with genomic positions and all info for each residues
    """

    logging.info(pdb_struct_id)

    # Ensembl mapping file
    e_mapping_df = (

        spark

        .createDataFrame(e_mapping_file[pdb_struct_id]['Ensembl'][gene_id]['mappings'])

        .select("chain_id", "genome_end", "genome_start", "end", "start")

        .withColumn("start_aut_resNb", f.col("start.author_residue_number"))
        .withColumn("end_aut_resNb", f.col("end.author_residue_number"))

        .drop("chains")
        .drop("end")
        .drop("start")

    )

    # Residues to obtain genomic position from
    res_info_df = spark.createDataFrame(residue_info)
    res_genomic_pos_df = (

        # Join residues with the mapping dataframe
        res_info_df
        .join(e_mapping_df,
              (res_info_df.protResNb >= e_mapping_df.start_aut_resNb) &
              (res_info_df.protResNb <= e_mapping_df.end_aut_resNb) &
              (res_info_df.chainId == e_mapping_df.chain_id)
              , "inner")

        .drop("chain_id")

        # Compute genomic positions of residues
        .withColumn("posRes1",
                        (((f.col("protResNb") - f.col("start_aut_resNb")) * 3) + f.col("genome_start")))
        .withColumn("posRes2",
                    (((f.col("protResNb") - f.col("start_aut_resNb")) * 3) + f.col("genome_start")) + 1)

        .withColumn("posRes3",
                    (((f.col("protResNb") - f.col("start_aut_resNb")) * 3) + f.col("genome_start")) + 2)

        # Remove entire row if residu position is already in the df
        .dropDuplicates(["posRes1", "chainId"])

        .groupBy("protResNb", "start_aut_resNb", "end_aut_resNb", "genome_start", "genome_end")

        .agg(f.collect_set(f.struct(f.col("protResNb"),
                                    f.col("intType"),
                                    f.col("chainId"),
                                    f.col("protResType"),
                                    f.col("posRes1"),
                                    f.col("posRes2"),
                                    f.col("posRes3")
                                    )
                           ).alias("ResInfos")
             )
    )

    # Store as a list for the apply (create a new column to the initial DF)
    all_residu__info_list = res_genomic_pos_df.select('ResInfos').collect()

    return all_residu__info_list


if __name__ == '__main__':

    program_description = '''
    Compute genomic position for each residues of each structure of the plip output.
    '''

    parser = argparse.ArgumentParser(add_help=True, description=program_description)

    parser.add_argument('-pi',
                        '--plip_input',
                        help='Path to the csv file with structures and drugs to compute PLIP interactions on.',
                        default=None,
                        metavar='csv_structure_drug_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-po',
                        '--plip_output',
                        help='Path to the output csv file with PLIP interactions computed for each structure-drug combinations.',
                        default=None,
                        metavar='csv_plip_output_file_path',
                        type=str,
                        required=True)

    parser.add_argument('-o',
                        '--output_folder',
                        help='Path to the output folder.',
                        default=None,
                        metavar='path_output_folder',
                        type=str,
                        required=True)

    parser.add_argument('-t',
                        '--test_set',
                        help='Add this argument to run the script on a small set (one structure)',
                        default=None,
                        action='store_true',
                        required=False)

    parser.add_argument('-l',
                        '--log_file',
                        help='File to save the logs to.',
                        default=None,
                        metavar='log_file',
                        type=str,
                        required=False)

    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, force=True)

    logging.info(program_description)

    main()
