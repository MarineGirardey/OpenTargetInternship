#!/bin/bash

# Download data
echo "Download all data"

# Drug OT
if [ ! -d "molecule" ]
then
    wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/parquet/molecule
fi

# Target OT
if [ ! -d "targets" ]
then
    wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/parquet/targets
fi

# Disease OT
if [ ! -d "diseases" ]
then
    wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/parquet/diseases
fi

# Target-disease evidence OT
if [ ! -d "evidence" ]
then
    wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/parquet/evidence
fi

# Inchikey crossref file
if [ ! -d "inchikey" ]
then
    mkdir inchikey
    cd inchikey
    wget https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/components_inchikeys.csv
    cd ..
fi

# Ensembl file to get target and chain id
if [ ! -d "ensembl" ]
then
    mkdir ensembl
    cd ensembl
    wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_ensembl.csv
    cd ..
fi

# Protein translation id mapping file to get the correct protein id
if [ ! -d "protein_mapping" ]
then
    mkdir protein_mapping
    cd protein_mapping
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
    cd ..
fi

# Mapping file for the correspondance between residue number and genomic position (generated by Daniel script)
if [ ! -d "residue_mapping" ]
then
    mkdir residue_mapping
    cd residue_mapping
    wget https://storage.googleapis.com/ot-team/marine/old_files/generated_mappings.tsv
    cd ..
fi

# Script 1
if [ ! -d "script_1_output" ]
then
    mkdir script_1_output
fi
time python OpenTargetInternship/scripts/get_struct_target_from_mol.py
  -m molecule
  -i inchikeys/components_inchikeys.csv
  -o script_1_output/structure_target_from_molecules/
  -e ensembl/pdb_chain_ensembl.csv
  -l script_1_output/structure_target_from_molecules/log.txt
  -f 
  2> errors.log

# Script 2
if [ ! -d "script_2_output" ]
then
    mkdir script_2_output
fi

if [ ! -d "script_2_output/structure_files" ]
then
    mkdir script_2_output/structure_files
fi

time python OpenTargetInternship/scripts/plip_interaction_mapping.py
  -i script_1_output/structure_drug_target_o.json/
  -o script_2_output/output_plip.json
  -l script_2_output/log.txt
  -f script_2_output/structure_files/ 
  2> errors.log

# Script 3
if [ ! -d "script_3_output" ]
then
    mkdir script_3_output
fi
time python OpenTargetInternship/scripts/residue_genomic_position.py 
  -o script_3_output/output_loc.json
  -i script_2_output/output_plip.json 
  -p protein_mapping/HUMAN_9606_idmapping.tsv
  -m residue_mapping/generated_mappings.tsv
  -l script_3_output/log.txt
  2> script_3_output/errors.log

# Script 4
if [ ! -d "script_4_output" ]
then
    mkdir script_4_output
fi
time python OpenTargetInternship/scripts/drug_disease_new_evidence.py
  -gl script_3_output/output_loc.json
  -e evidence/
  -m molecule
  -i inchikey/components_inchikeys.csv
  -d diseases/
  -t targets/
  -o script_4_output/output_evidence_filtered.csv
  2> script_4_output/errors.log