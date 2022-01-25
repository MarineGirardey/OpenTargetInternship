import pyspark.sql.functions as F
from pyspark.sql import SparkSession
import os
import requests
import timeit
global spark

'''
This script retrieves the structures in complex with the molecules from parquet.
Output example:

ChEmbl ID | PDB_ID | struct
---------------------------

'''


def convert_molecule_file_to_dataframe(path_molecule):
    """
    Insert all molecules form parquet files in a unique df
    --------------
    Keyword arguments
        path_molecule: path to molecule parquet files
    --------------
    Return
        molecule : a df containing molecules
    """

    # it for the first molecule insertion
    it = 0

    for molecule_filename in os.listdir(path_molecule):
        # Store id column of molecule file in a temporary df
        temp_df = (
            spark.read
            .parquet(path_molecule + molecule_filename, header=True)
            .select(F.col("id"))
        )

        # Because union doesn't work for empty df
        if it == 0:
            molecule_df = temp_df

        # Add rows from temp_df to the all files df "molecule_df"
        else:
            molecule_df = molecule_df.union(temp_df)

        it += 1

    return molecule_df


def join_molecule_and_unichem_df(molecule_df, unichem_df):
    """
    Join df with all molecules and df with UniChem ID in one df
    --------------
    Keyword arguments
        molecule_df: df with all molecules from parquet files
        unichem_df: df with molecules ID from UniChem
    --------------
    Return
        molecule_in_unichem_df : a df with only molecules from parquet retrieved in UniChem (with UniChem ID)
    """

    molecule_in_unichem_df = (molecule_df
                              .join(unichem_df, unichem_df["From src:'1'"] == molecule_df["id"])
                              .withColumnRenamed("id", "MOLECULE_ID")
                              .withColumnRenamed("From src:'1'", "CHEMBL_ID")
                              .withColumnRenamed("To src:'3'", "PDB_ID")
                              )
    molecule_in_unichem_df = molecule_in_unichem_df.drop('MOLECULE_ID')

    return molecule_in_unichem_df


def get_structure(pdb_id):
    """
    Scrap the list of structures in complex with each molecule from PDBe API
    Function to apply to a df column
    --------------
    Keyword arguments
        pdb_id: id of the molecule for whom we want to know the structures associated
    --------------
    Return
        data[pdb_id] : data is a dictionary, we want the value which is the list of structures
    """

    url = f'https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/{pdb_id}'
    response = requests.get(url)
    try:
        data = response.json()
        return data[pdb_id]
    except:
        return None


def apply_get_structure_function_on_df(unichem_molecule_joined_df):
    """
    Apply on a df column the get_structure function to associate for each PDB ID the structures in complex with
    --------------
    Keyword arguments
        unichem_molecule_joined_df: df with UniChem ID and PDB ID
    --------------
    Return
        unichem_molecule_struct_pd_df: pandas df with UniChem ID, PDB ID and list of structures
    """

    unichem_molecule_struct_pd_df = (unichem_molecule_joined_df
                                     .toPandas()
                                     .assign(struct=lambda x: x['PDB_ID'].apply(get_structure)))
    return unichem_molecule_struct_pd_df


def all_statistics(unichem_molecule_df, unichem_molecule_struct_pd_df, unichem_molecule_struct_spark_df, molecule_df, unichem_df):
    """
    Compute statistics for our knowledge
    --------------
    Keyword arguments
        unichem_molecule_df: df with UniChem ID and molecule ID
        unichem_molecule_struct_pd_df: pandas df with UniChem ID, molecule ID and structures
        unichem_molecule_struct_spark_df: same but spark df
        molecule_df: df only with molecule ID
        unichem_df: df only with UniChem ID
    --------------
    Return
        total_nb_unichem: number of molecule in UniChem BDD
        total_nb_molecule: number of molecule in all parquet files
        nb_molecule_in_unichem: number of our molecule found in UniChem
        percentage_molecule_in_unichem: percentage of molecules in UniChem
        total_nb_struct: total number of structure
        nb_molecule_with_struct: number of molecule with structures
    """

    total_nb_unichem = unichem_df.count()
    total_nb_molecule = molecule_df.count()
    nb_molecule_in_unichem = unichem_molecule_df.count()
    percentage_molecule_in_unichem = (nb_molecule_in_unichem / total_nb_molecule) * 100
    total_nb_struct = len(unichem_molecule_struct_pd_df.explode('struct'))
    nb_molecule_with_struct = unichem_molecule_struct_spark_df.count()

    return total_nb_unichem, total_nb_molecule, nb_molecule_in_unichem, percentage_molecule_in_unichem, total_nb_struct, nb_molecule_with_struct


def main():
    start1 = timeit.default_timer()
    spark = SparkSession.builder.getOrCreate()

    path_molecule_files = '/molecule/'
    path_unichem_file = '/id_files/src1src3.txt'

    unichem_df = spark.read.csv(path_unichem_file, sep=r'\t', header=True)

    molecule_df = convert_molecule_file_to_dataframe(path_molecule_files)
    unichem_molecule_joined_df = join_molecule_and_unichem_df(molecule_df, unichem_df)

    start2 = timeit.default_timer()
    unichem_molecule_struct_pd_df = apply_get_structure_function_on_df(unichem_molecule_joined_df)
    unichem_molecule_struct_spark_df = spark.createDataFrame(unichem_molecule_struct_pd_df)
    stop2 = timeit.default_timer()

    unichem_molecule_struct_pd_df.to_csv("structure_of_molecules.csv", index=False, header=True)

    statistics = all_statistics(unichem_molecule_joined_df,
                                unichem_molecule_struct_pd_df,
                                unichem_molecule_struct_spark_df,
                                molecule_df,
                                unichem_df)

    stop1 = timeit.default_timer()

    return statistics, start1, stop1, start2, stop2


if __name__ == '__main__':
    main = main()
    stats = main[0]
    start1 = main[1]
    stop1 = main[2]
    start2 = main[3]
    stop2 = main[4]

    print('Total time running of the script: ', round(stop1 - start1, 2))
    print('With ', round(stop2 - start2, 2), 'of time for scrap the structure for each molecule on the API')

    print("Counting: %s molecules in UniChem in total" % stats[0])
    print("Counting: %s molecules from parquet in total" % stats[1])
    print("Counting: %s of our molecules in UniChem" % stats[2])
    print("%s of molecules are in UniChem" % (str(round(stats[3], 2)) + '%'))
    print("Counting: %s structures int total" % stats[4])
    print("Counting: %s molecules with structure(s)" % stats[5])
