#!/usr/bin/env bash
set -Eeu

# Script based on AQUAMIS Setup Script 
# https://gitlab.com/bfr_bioinformatics/AQUAMIS/-/blob/master/scripts/aquamis_setup.sh
# BSD 3-Clause License

# Copyright (c) 2019, Carlus Deneke, Holger Brendebach and Simon H. Tausch
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
  # list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
  # this list of conditions and the following disclaimer in the documentation
  # and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
  # contributors may be used to endorse or promote products derived from
  # this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# retrieve databases and extend envs to use in aquamis

msg() {
  echo >&2 -e "${1-}"
}

die() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  exit "$code"
}

arg_databases=1
arg_busco=1
set -x

# Retrieve script paths
SCRIPT_PATH=${CONDA_PREFIX}/opt/aquamis/scripts/
INSTALL_PATH=$(dirname $SCRIPT_PATH)

# Fancy echo
bold=$(tput bold);normal=$(tput sgr0);red=$(tput setaf 1);green=$(tput setaf 2);yellow=$(tput setaf 3);blue=$(tput setaf 4);magenta=$(tput setaf 5);cyan=$(tput setaf 6);grey=$(tput setaf 7);inverted=$(tput rev)


# Download Function
download_file() {
  download_success=0
  download_hash=''
  local target=$1
  local source=$2
  if [[ -f $target ]] && [[ $arg_force == 0 ]]; then
    echo "The file $target already exists. Skipping download."
  else
    echo "Downloading $source to $target"
    wget --ca-certificate=${INSTALL_PATH}/resources/seafile-bfr-berlin_certificate.pem --output-document $target $source
    [[ $? -eq 0 ]] && [[ -s $target ]] && download_hash=$(openssl dgst -r -sha256 $target) && download_success=1
  fi
}

# Download Databases
complete_busco() {
  
  echo "Downloading Augustus and BUSCO databases to ${INSTALL_PATH}/download and extracting them in ${CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/"

  # Create Subdirectories
  [[ ! -d "${INSTALL_PATH}/download" ]] && mkdir -p ${INSTALL_PATH}/download
  cd ${INSTALL_PATH}/download

  # Download files
  [[ ! -d "${CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/" ]] && die "ERROR: The Conda environment installation is incomplete. QUAST was not installed. Aborting..."
  download_file "augustus.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/augustus.tar.gz" # 139MB or https://databay.bfrlab.de/f/fa13c9eb2625477eb729/?dl=1
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f augustus.tar.gz -C ${CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/

  download_file "bacteria.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/bacteria.tar.gz" # 8.8MB or https://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f bacteria.tar.gz -C ${CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/busco/

  echo "${green}BUSCO installation completed.${normal}"
}

download_databases() {
  echo "Downloading databases to ${INSTALL_PATH}/download and extracting them in their respective folder under ${INSTALL_PATH}/reference_db"

  # Create Subdirectories
  [[ ! -d "${INSTALL_PATH}/download" ]] && mkdir -p ${INSTALL_PATH}/download
  mkdir -p ${INSTALL_PATH}/reference_db/mash/genomes
  mkdir -p ${INSTALL_PATH}/reference_db/taxonkit
  cd ${INSTALL_PATH}/download

  # Download files
  download_file "confindr_db.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/confindr_db.tar.gz" # 153MB for the old db incl. some rMLST db's e.g. Enterococcus and Brucella, use https://seafile.bfr.berlin/f/47cae689eda7440c83bb/?dl=1
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f confindr_db.tar.gz -C ${INSTALL_PATH}/reference_db/

  download_file "mashDB.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/mashDB.tar.gz" # 216MB or https://databay.bfrlab.de/f/34b8c88945a8439dac64/?dl=1
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f mashDB.tar.gz -C ${INSTALL_PATH}/reference_db/mash/

  download_file "minikraken2.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/minikraken2.tar.gz" # 6.0GB or ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz or https://genome-idx.s3.amazonaws.com/kraken/minikraken_8GB_202003.tgz
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f minikraken2.tar.gz -C ${INSTALL_PATH}/reference_db/

  download_file "taxdump.tar.gz" "https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/taxdump.tar.gz" #  54MB or ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
  [[ -n $download_hash ]] && echo "$download_hash" >> reference_db.sha256
  [[ "$download_success" == 1 ]] && tar -xzv -f taxdump.tar.gz -C ${INSTALL_PATH}/reference_db/taxonkit/

  echo "${green}Download finished.${normal}"
}

## Execute Option Modules
[[ "$arg_busco" == 1 ]] && complete_busco
[[ "$arg_databases" == 1 ]] && download_databases
