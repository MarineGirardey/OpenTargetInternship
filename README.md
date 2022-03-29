# OpenTargetInternship
6 months internship at the EMBL-EBI in the data team of the OpenTarget platform

## Create Google Compute Engine instance
```bash
# Set parameters.
export INSTANCE_NAME=plip-interaction
export INSTANCE_ZONE=europe-west1-d
export INSTANCE_TYPE=n1-highcpu-64

# Create the instance and SSH.
gcloud compute instances create \
  ${INSTANCE_NAME} \
  --project=open-targets-eu-dev \
  --zone=${INSTANCE_ZONE} \
  --machine-type=${INSTANCE_TYPE} \
  --service-account=426265110888-compute@developer.gserviceaccount.com \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --create-disk=auto-delete=yes,boot=yes,device-name=${INSTANCE_NAME},image=projects/ubuntu-os-cloud/global/images/ubuntu-2004-focal-v20210927,mode=rw,size=2000,type=projects/open-targets-eu-dev/zones/europe-west1-d/diskTypes/pd-balanced

# SSH command may take a while to work while the instance is provisioned and configured.
gcloud compute ssh --zone ${INSTANCE_ZONE} ${INSTANCE_NAME}

# Use screen to avoid losing output when connection is lost. On reconnect, the session can be restored with calling `screen -d -r`.
screen
```

## Configure environment and install dependencies (first time)
```bash
sudo apt update
# python3-openbabel has to be installed globally because of a number of errors in the current PIP packaging
# https://github.com/openbabel/openbabel/issues/2408
sudo apt install -y \
  openbabel \
  python3-openbabel \
  python3-pip \
  python3-testresources \
  python3.8-venv

git clone -q https://github.com/MarineGirardey/OpenTargetInternship
cd OpenTargetInternship
python3 -m venv --system-site-packages venv
source venv/bin/activate
pip install --quiet --upgrade pip setuptools
# Manually mark openbabel as installed, because we did it previously and we don't want pip to try to install the broken PyPi version
# https://stackoverflow.com/questions/39403002/manually-set-package-as-installed-in-python-pip
touch venv/lib/python3.8/site-packages/openbabel-3.0.0-py3.8.egg-info
pip install --quiet --upgrade \
  dask \
  distributed \
  matplotlib \
  numpy \
  pandarallel \
  pandas \
  plip \
  pyspark \
  requests \
  git+https://github.com/PDBeurope/arpeggio
```

## Commands to reconnect to the machine and/or reactivate the environment
* Reconnect: `gcloud compute ssh --zone ${INSTANCE_ZONE} ${INSTANCE_NAME}`
* Restore previously created screen session: `screen -d -r`
* Reactivate the environment: `cd ~/OpenTargetsInternship && source venv/bin/activate`

## Run the analysis
```bash
cd scripts
mkdir -p pdb
time python script_plip_interaction_mapping.py \
  --input_file structure_for_plip.csv \
  --output_file output.csv \
  --log_file log.txt \
  --pdb_folder pdb
```
