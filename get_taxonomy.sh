#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a get_taxonomy.`date +%Y-%m-%d`.log)
exec 2> >(tee -a get_taxonomy.`date +%Y-%m-%d`.log >&2)

##### INFO

# get_taxonomy.sh

# for generating the taxonomy file for protax (using Integrated Taxonomic Information System)

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash get_taxonomy.sh taxonName taxonRank screenforbio
# where:
# taxonName is the scientific name of the target taxon e.g. Tetrapoda
# taxonRank is the classification rank of the target taxon e.g. superclass
# screenforbio is the path to the screenforbio-mbc directory

# note: only ITIS classification for kingdom Animalia supported at present

##### WELCOME
if [ "$#" -ne 3 ]
then
	echo ""
	echo "You are trying to use get_taxonomy.sh but have not provided enough information."
	echo "Please check the following:"
	echo ""
	echo "Usage: bash get_taxonomy.sh taxonName taxonRank screenforbio"
	echo "Where:"
	echo "taxonName is the scientific name of the target taxon e.g. Tetrapoda"
	echo "taxonRank is the classification rank of the target taxon e.g. superclass"
	echo "screenforbio is the path to the screenforbio-mbc folder"
	echo ""
	exit 1
else
	echo ""
	echo "Welcome to get_taxonomy.sh: fetching ITIS taxonomy for target taxon"
	echo ""
	echo "The time now is $(date)."
	echo ""
	echo "Looking up taxonomy for ${2} ${1}"
	echo ""
fi

##### PARAMETERS
TAXON=${1}
RANK=${2}
SCRIPTS=${3}

##### MAIN
echo "Running get_taxonomy.R now - this may take some time..."
echo ""
start=`date +%s`
Rscript --vanilla ${SCRIPTS}/get_taxonomy.R ${TAXON} ${RANK}
echo ""
echo "get_taxonomy.R finished at $(date)"
end=`date +%s`
runtime=$((end-start))
echo ""

##### END
echo "Fetching taxonomy of ${TAXON} took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
echo ""
echo "Next step: Make curated non-redundant reference database with get_sequences.sh"
echo ""
echo "Have a nice day :-)"
echo ""
