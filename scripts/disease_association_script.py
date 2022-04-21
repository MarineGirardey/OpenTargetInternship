import pyspark.sql.functions as f
from pyspark.sql import SparkSession
spark = SparkSession.builder.getOrCreate()
from pyspark.sql.types import StructField, StructType, IntegerType, StringType


evidence = (
    spark.read.parquet("/home/marinegirardey/evidence")
    .select("datasourceId", "targetId", "variantId", "diseaseId", "diseaseFromSource")
    .filter(f.col("variantId").isNotNull())
    .withColumn('chr', f.split(f.col('variantId'), '_').getItem(0))
    .withColumn('genomicLocation', f.split(f.col('variantId'), '_').getItem(1))
)

gen_location = (
    spark.read.json("residue_gen_pos_output/residue_genomic_position_2.json")
)

gen_location = (
    gen_location
    .withColumn("pdbCompound", gen_location["resInfos.compound"])
    .withColumn("resNb", gen_location["resInfos.res_nb"])
    .withColumn("chain", gen_location["resInfos.chain"])
    .withColumn("resType", gen_location["resInfos.res_type"])
    .withColumn("interType", gen_location["resInfos.inter_type"])
    .withColumn("chr", gen_location["resInfos.chromosome"])
    .withColumn("genLocation_1", gen_location["resInfos.genLocation.res_pos_1"])
    .withColumn("genLocation_2", gen_location["resInfos.genLocation.res_pos_2"])
    .withColumn("genLocation_3", gen_location["resInfos.genLocation.res_pos_3"])
    .drop("resInfos")
)

res_with_disease = (
    gen_location.join(
        evidence, 
        (gen_location.chr == evidence.chr) &
        (
            (gen_location.genLocation_1 == evidence.genomicLocation) |
            (gen_location.genLocation_2 == evidence.genomicLocation) |
            (gen_location.genLocation_3 == evidence.genomicLocation)
        )
    )
)

print(res_with_disease.count())
