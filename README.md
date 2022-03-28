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
  matplotlib \
  numpy \
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

```