#!/usr/bin/env bash
set -Eeu

# simply copies the prodigal folder from chewiesnake env to 
# Repo basedir

PRODIGAL=${CONDA_PREFIX}/opt/chewieSnake/chewBBACA/CHEWBBACA/prodigal_training_files

SCRIPT_DIR=$( dirname -- "$0"; )
WF=$(dirname ${SCRIPT_DIR})
REPO=$({dirname $WF})

cp -r ${PRODIGAL} ${REPO}