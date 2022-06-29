import pyspark.sql.functions as f
from pyspark.sql import SparkSession
from pyspark.sql.types import StructField, StructType, IntegerType, StringType
import argparse
import logging
import pandas as pd
import sys
# Get residue nb and type
import io
import requests


def main():

    spark = (
        SparkSession.builder
        .master('local[*]')
        .config("spark.driver.memory", "15g")
        .getOrCreate()
    )

    # Open genomic location output file
    res_gen_loc = (
        spark.read.json(args.gen_loc)
        .withColumnRenamed("pdb_structure_id", "pdbStructureId")
        .withColumnRenamed("prot_chain_id", "chainId")
        .withColumnRenamed("compound_id", "pdbCompoundId")
    )

    # Pivot expression to reduce pyspark operation
    unpivot_exp = '''stack(3, 'pos1', pos1, 'pos2', pos2, 'pos3', pos3) as (genLocation_label, genLocation_val)'''

    # Select column and apply the unpivot expression
    gen_location_unpivot = (
        res_gen_loc
        .select('geneId', 'accession', 'pdbStructureId',
                'pdbCompoundId', 'prot_residue_number',
                'chainId', 'prot_residue_type', 'chr',
                f.expr(unpivot_exp))
    )

    # OT platform target-disease dataset (variant id - disease association)
    evidence = (
        spark.read.parquet(args.evidence)
        .filter(f.col("variantId").isNotNull())
        .withColumn('chromosome', f.split(f.col('variantId'), '_').getItem(0))
        .withColumn('genomicLocation', f.split(f.col('variantId'), '_').getItem(1))
        .groupBy('chromosome', 'genomicLocation')
        .agg(
            f.collect_set(f.struct(
                f.col('variantId'),
                f.col('diseaseId'),
                f.col('diseaseFromSource'),
                f.col('datasourceId'))).alias('evidenceInfo')
        )
    )

    # Join on chromosome and residue genomic location with variant id to search residues associated with a disease
    disease_association = (
        gen_location_unpivot
        .join(
            evidence,
            (
                gen_location_unpivot.chr == evidence.chromosome
            ) & (
                gen_location_unpivot.genLocation_val == evidence.genomicLocation
            ),
            how='left'
        )
    )

    # Rearrange the dataset
    associations = (
        disease_association
        .select('*', f.explode('evidenceInfo').alias('evidence'))
        .select('*', f.col('evidence.*'))
        .select(
            f.col('geneId').alias('targetId'),
            f.col('pdbCompoundId').alias('pdbCompound'),
            'diseaseId',
            'variantId'
        )
        .distinct()
        .persist()
    )

    # Drug dataset with inchikey id
    molecule = (
        spark.read
        .parquet(args.molecule)
        .select(
            f.col('inchiKey').alias('inchikey'),
            f.col('name').alias('drugName')
        )
        .persist()
    )

    # Inchikey - pdb id correspondance file
    inchikey = (
        spark.read
        .csv(args.inchikey, sep=',', header=True, comment='#')
        .select(
            f.col('InChIKey').alias('inchikey'),
            f.col('CCD_ID').alias('pdbCompound')
        )
        .persist()
    )

    # Drug dataset with compound id
    molecules_inchikey_join = (
        molecule
        .join(inchikey, on='inchikey')
        .drop("inchikey")
        .persist()
    )

    # Target dataset
    targets = (
        spark.read
        .parquet(args.target)
        .select(
            f.col('id').alias('targetId'),
            f.col('approvedSymbol').alias('symbol'),
            f.col('approvedName').alias('targetName')
        )
        .persist()
    )

    # Disease dataset
    diseases = (
        spark.read
        .parquet(args.disease)
        .select(
            f.col('id').alias('diseaseId'),
            f.col('name').alias('diseaseName')
        )
        .persist()
    )

    # Join data to obtain the target, drug, disease name info
    mapped_associations = (
        associations
        .groupBy('targetId', 'pdbCompound', 'diseaseId')
        .agg(f.collect_set('variantId').alias('variantIds'))

        # Joining with disease name
        .join(diseases, on='diseaseId', how='left')

        # Joining with gene name
        .join(targets, on='targetId', how='left')

        # Joining with drug name
        .join(molecules_inchikey_join, on='pdbCompound')

        .persist()
    )

    # Filter out unwanted compounds not removed with biolip and reorganise
    disease_evidence = (
        mapped_associations
        .filter(
            (f.col('pdbCompound') != 'ATP')
            & (f.col('pdbCompound') != 'ZN')
        )
        .drop('targetId', 'diseaseId')
        .withColumnRenamed('name', 'targetName')
        .na.drop("any")
    )

    # TODO: Avoid this
    # Re open molecule
    molecule = (
        spark.read.parquet(args.molecule)
        .select(
            f.col('id'),
            f.col('maximumClinicalTrialPhase'),
            f.col('inchiKey').alias('inchikey'),
            f.col('name').alias('drugName'),
            'linkedTargets', 'linkedDiseases'
        )
    )

    # Join molecule and target dataset (For linkedTargetName)
    drugs_w_linked_target_names = (
        molecule
        .filter(f.col('linkedTargets').isNotNull())
        .select(
            f.col('id').alias('drugId'),
            f.col('maximumClinicalTrialPhase'),
            f.explode(f.col('linkedTargets.rows')).alias('targetId')
        )

        .join(targets, on='targetId', how='left')
        .groupby('drugId')
        .agg(
            f.collect_set('targetName').alias('linkedTargetName')
        )
        .persist()
    )

    # Join molecule and disease (For linkedDiseasesName)
    drugs_w_linked_disease_names = (
        molecule
        # .filter(f.col('linkedDiseases').isNotNull())
        .select(
            f.col('id').alias('drugId'),
            f.col('maximumClinicalTrialPhase'),
            f.explode(f.col('linkedDiseases.rows')).alias('diseaseId')
        )
        .join(diseases, on='diseaseId', how='left')
        .groupby('drugId')
        .agg(
            f.collect_set('diseaseName').alias('linkedDiseaseName')
        )
        .persist()
    )

    # Join molecule and molecule+target and molecule+disease
    resolved_molecules = (
        molecule
        .select(
            f.col('id').alias('drugId'),
            f.col('maximumClinicalTrialPhase'),
            f.col('drugName').alias('drugName'),
        )
        .join(drugs_w_linked_disease_names, on='drugId', how='left')
        .join(drugs_w_linked_target_names, on='drugId', how='left')
    )

    # Join previous dataset with the infos with the evidence dataset
    joined_mapped_linked = (
        disease_evidence
        .select(
            f.col('drugName'),
            f.col('diseaseName').alias('evidenceDiseaseName'),
            f.col('targetName').alias('interactingTargetName'),
            f.col('variantIds')
        )
        .join(resolved_molecules, on='drugName', how='left')
        .persist()
    )

    # Collecting list of variants
    v_list = (
        joined_mapped_linked
        .filter(f.col('linkedDiseaseName').isNotNull())
        .select(f.explode('variantIds').alias('variantId'))
        .distinct()
        .collect()
    )
    variant_list = [variant.variantId.replace('_', ' ') for variant in v_list]

    # Get residue info from ProtVar
    protVarMappings_pdf = get_variant_info_protvar(variant_list)

    # Dict with variant and residue nb and type in the good format
    protVar_map = (
        protVarMappings_pdf
        .assign(
            variantId=lambda df: df.User_input.str.replace(' ', '_'),
            mutation=lambda df: df.apply(get_mutation, axis=1)
        )
        [['variantId', 'mutation']]
        .groupby('variantId')
        .agg({
            'mutation': lambda s: ' '.join(s.loc[s.notna()].to_list())
        })
        .mutation
        .to_dict()
    )

    map_mutations = f.udf(
        lambda variantIds: ' '.join([protVar_map[variant] for variant in variantIds if variant in protVar_map]),
        StringType()
    )

    # Dataset with the residue number and type info for each mutation
    w_mutations = (
        joined_mapped_linked
        .filter(f.col('linkedDiseaseName').isNotNull())
        .withColumn('mutations', map_mutations(f.col('variantIds')))
        .persist()
    )

    # Filter on the approval phase, remove cancer cases and reorganise
    w_mutations_filtered = (
        w_mutations

        .filter(f.col('maximumClinicalTrialPhase') >= 4)
        .select(f.explode(f.col('linkedDiseaseName')).alias('linkedDiseaseNameExploded'), '*')
        .filter(~f.col('linkedDiseaseNameExploded').contains('cancer'))
        .filter(~f.col('evidenceDiseaseName').contains('cancer'))

        .drop('linkedDiseaseName')

        .groupby(
            f.col('drugName'),
            f.col('evidenceDiseaseName'),
            f.col('interactingTargetName'),
            f.col('drugId'),
            f.col('maximumClinicalTrialPhase'),
            f.col('linkedTargetName'),
            f.col('mutations')
        )
        .agg(
            f.collect_set(f.col('linkedDiseaseNameExploded')).alias('linkedDiseaseName')
        )
    )

    # Reorganise and converst to pandas
    w_mutations_filt_agg = (
        w_mutations_filtered
        .groupBy('drugName', 'interactingTargetName')
        .agg(
            f.collect_set('evidenceDiseaseName').alias('evidenceDiseaseNames'),
            f.first('drugId'),
            f.first('linkedDiseaseName'),
            f.first('linkedTargetName'),
            f.collect_set('mutations')
        )
        .toPandas()
    )

    # TODO: Add this line (need to be corrected)
    # interesting_joined_mapped_linked = (
    #     w_mutations_filt_agg
    #     .toPandas()
    #     .assign(
    #         new_values=lambda df: df.apply(
    #             lambda row: [
    #                 a for a in row['first(linkedTargetName)'] if a not in row['interactingTargetName']
    #             ], axis=1)
    #     )
    # )

    print(w_mutations_filt_agg)

    # Save result in a csv (because small dataset and to read it easely with Excel)
    w_mutations_filt_agg.to_csv(args.output)

def get_mutation(row: tuple) -> str:
    """Get only number and type of the residue to create a new column with this info
    Args:
        row: row of the dataframe
    Returns:
        A pandas DataFrame all info on the mutations from ProtVar
    """

    try:
        position = int(row['Amino_acid_position'])
        return row['Amino_acid_change'].replace('/', str(position))

    except:
        return None

def get_variant_info_protvar(variant_id: list) -> pd.DataFrame:
    """Get mutation information from ProtVar to obtain the number and AminoAcid type of the mutation
    Args:
        variant_id: list with the variant id of the new drug-disease aqssociation found
    Returns:
        A pandas DataFrame all info on the mutations from ProtVar
    """

    # ProtVar API
    URL = 'https://www.ebi.ac.uk/ProtVar/api/download/stream'

    headers = {
        'accept': '*/*',
    }

    params = {
        'function': 'false',
        'population': 'false',
        'structure': 'false',
    }

    # Request to get API response
    response = requests.post(URL, params=params, headers=headers, json=variant_id)

    # Output
    return (
        pd.read_csv(io.StringIO(response.content.decode('utf-8')))
        [[
            'User_input',
            'Gene',
            'Codon_change',
            'Amino_acid_change',
            'Protein_name',
            'Amino_acid_position',
            'Consequences'
        ]]
    )


if __name__ == '__main__':

    global spark

    program_description = '''
    This script join and filter data in order to investigate potential new Drug-Diseases associations.
    For all OT files : https://platform.opentargets.org/downloads
    inchikey file : https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/components_inchikeys.csv
    '''

    parser = argparse.ArgumentParser(add_help=True, description=program_description)

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(
        description='This script maps plip output to Ensembl translations.')

    parser.add_argument('--gen_loc', '-gl', type=str,
                        help='Input file with the genomic locations.', required=True)

    parser.add_argument('--evidence', '-e', type=str,
                        help='target-diseases evidence from OT platform.', required=True)

    parser.add_argument('--molecule', '-m', type=str,
                        help='molecule from OT platform.', required=False)

    parser.add_argument('--target', '-t', type=str,
                        help='target from OT platform.', required=True)

    parser.add_argument('--inchikey', '-i', type=str,
                        help='inchikey cross ref for compound ids.', required=True)

    parser.add_argument('--disease', '-d', type=str,
                        help='target from OT platform.', required=True)

    parser.add_argument('--output', '-o', type=str,
                        help='output file with new potential target_disease evidence.', required=True)

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        # filename=args.log_file,
        force=True)

    logging.StreamHandler(sys.stderr)

    logging.info(program_description)

    main()
