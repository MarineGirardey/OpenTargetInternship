# Enrichment of the Open Targets Platform with structural annotations

Marine Girardey - marine.girardey@gmail.com

<a href="https://www.linkedin.com/in/marine-girardey/">
  <img src="https://img.shields.io/badge/linkedin-%230077B5.svg?&style=for-the-badge&logo=linkedin&logoColor=white" />
</a>&nbsp;&nbsp;

[DLAD Master](https://formations.univ-amu.fr/fr/master/5SBI/PRSBI5AB) - [Aix-Marseille University](https://www.univ-amu.fr/en)

6 months master internship at the [EMBL-EBI](https://www.ebi.ac.uk/), Hinxton in the data team of [Open Targets](https://www.opentargets.org/)

<p align='center'>
  <img width="387" alt="ot" src="https://user-images.githubusercontent.com/77064839/176146051-e7d298d7-7863-4a12-978f-6514f6cd8eb3.png">
  <img width="335" alt="ot" src="https://user-images.githubusercontent.com/77064839/176423103-8755774d-af4c-4997-8804-5ae0854d1096.png">
</p>

## Context

The Open Targets Platform provide annotations for target prioritization but none of the annotation comes from structural information.

The goal of this project is to enrich the platform with structural information about the drug-target complex to provide an interactive display 
of the 3D complex on the platform and create a new dataset from a new structure-based new association investigation.
For more details, read my [internship report]().

### How to create Google Compute Engine instance?
```bash
# Set parameters.
export INSTANCE_NAME=mgirardey_project
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

### How to configure the environment needed?
```bash
sudo apt update
# python3-openbabel has to be installed globally because of a number of errors in the current PIP packaging
# https://github.com/openbabel/openbabel/issues/2408
sudo apt install -y \
  openbabel \
  openjdk-11-jre-headless \
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
  distributed \
  matplotlib \
  numpy \
  pandarallel \
  pandas \
  plip \
  pyspark \
  requests \
```

### Commands to reconnect to the machine and/or reactivate the environment
* Reconnect: `gcloud compute ssh --zone ${INSTANCE_ZONE} ${INSTANCE_NAME}`
* Restore previously created screen session: `screen -d -r`
* Reactivate the environment: `cd ~/OpenTargetInternship && source venv/bin/activate`


### How to run the full code?
### How to run each script independently?

# PLIP


# Get Location

## Run the analysis
```bash
cd scripts
mkdir -p pdb
time python script_plip_interaction_mapping.py \
  --input_file structure_for_plip_human_structures.csv \
  --output_file output.csv \
  --log_file log.txt \
  --pdb_folder pdb
mkdir -p residue_gen_pos_output
time python residue_genomic_position_script.py \
  --plip_input gene_mapped_structures.json \
  --plip_output output.csv \
  --output_folder residue_gen_pos_output \
  --log_file genomic_position_log.txt
```


python scripts/residue_genomic_position.py -o output_file_location/output_loc.json -i output_file_plip/output_plip.json -p scripts/files_to_merge_genomic_loc/HUMAN_9606_idmapping.tsv -m scripts/files_to_merge_genomic_loc/generated_mappings.tsv -i output_file_plip/output_plip.json -l output_file_location/log.txt
