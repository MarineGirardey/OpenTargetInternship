# OpenTargetInternship
6 months internship at the EMBL-EBI in the data team of the OpenTarget platform

## Configure environment and install dependencies (first time)
```bash
sudo apt update
sudo apt install -y \
  openbabel \
  python3-openbabel \
  python3.8-venv
git clone -q https://github.com/MarineGirardey/OpenTargetInternship
cd OpenTargetInternship
python3 -m venv venv
source venv/bin/activate
pip install --quiet --upgrade pip setuptools
pip install --quiet --upgrade \
  dask \
  distributed \
  matplotlib \
  numpy \
  pandarallel \
  pandas \
  pyspark \
  requests \
  git+https://github.com/PDBeurope/arpeggio \
  git+https://github.com/pharmai/plip
```

## Reactivate the prepared environment
```bash
source venv/bin/activate
```

## Run the analysis
```bash
cd scripts
mkdir pdb
time python script_plip_interaction_mapping.py \
  --input_file structure_for_plip_small_set.csv \
  --output_file structure_for_plip_small_set_output.csv\
  --pdb_folder pdb \
  --nb_partitions 600
```
