#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a train_weighted_protax.`date +%Y-%m-%d`.log)
exec 2> >(tee -a train_weighted_protax.`date +%Y-%m-%d`.log >&2)

##### INFO

# train_weighted_protax.sh

# for selecting training data, running required LAST searches and protax model parameterisation

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash train_weighted_protax.sh splist taxonomy screenforbio
# where:
# splist is a list of expected species to use in weighing in the format Genus,species (e.g. Homo,sapiens)
# taxonomy is the final protax-formatted taxonomy file from get_sequences.sh (e.g. Tetrapoda.final_protax_taxonomy.txt)
# screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)

# note: will take the taxon from the protax taxonomy file name
# note: assumes curated database FASTA files are in current directory and labelled with format taxon.final_database.locus.fa

##### WELCOME
if [ "$#" -ne 3 ]
then
  echo ""
	echo "You are trying to use train_weighted_protax.sh but have not provided enough information."
	echo "Please check the following:"
	echo ""
	echo "usage: bash train_weighted_protax.sh splist taxon screenforbio"
	echo "where:"
  echo "splist is a list of expected species to use in weighing in the format Genus,species (e.g. Homo,sapiens)"
  echo "taxonomy is the final protax-formatted taxonomy file from get_sequences.sh (e.g. Tetrapoda.final_protax_taxonomy.txt)"
  echo "screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)"
  exit 1
else
  echo "Welcome to train_weighted_protax.sh: preparing for and parameterising weighted PROTAX models"
  echo ""
  echo "The time now is $(date)."
  echo "Information from you:"
  echo "Species for weighting are in ${1}"
  echo "Taxonomy file is in ${2}"
  echo "Required Perl scripts are in ${3}/protaxscripts"
  echo ""
fi

##### PARAMETERS
SPLIST=${1}
TAXONOMY=${2}
SCRIPTS=${3}
DATE=`date +%Y-%m-%d`
taxon=$(basename ${TAXONOMY} | cut -f1 -d".")

##### MAIN
start=`date +%s`
# get file names
if ls ${taxon}.final_database.*.fa 1> /dev/null 2>&1
then
FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "${taxon}.final_database.*.fa" | sed 's/\.\///g'))
else
  echo "No matching FASTA databases found."
  echo ""
  exit 1
fi
echo "Step 1: make taxonomy priors"
for file in ${FILES[@]}
do
  locus=$(basename $file ".fa" | cut -f3 -d".")
  echo "Working on ${locus}"
  mkdir ./w_model_${locus}
  perl ${SCRIPTS}/protaxscripts/maketaxonomy.pl ${TAXONOMY} > ./w_model_${locus}/taxonomy
  perl ${SCRIPTS}/protaxscripts/taxonomy_priors.pl ./w_model_${locus}/taxonomy > ./w_model_${locus}/tax4
  grep -f ${SPLIST} ./w_model_${locus}/tax4 > ./w_model_${locus}/expected_sp
  perl ${SCRIPTS}/protaxscripts/setpriors.pl 0.9 ./w_model_${locus}/expected_sp ./w_model_${locus}/tax4 > ./w_model_${locus}/wtax4
  rm ./w_model_${locus}/tax4
  rm ./w_model_${locus}/expected_sp
  perl ${SCRIPTS}/protaxscripts/thintaxonomy.pl 1 ./w_model_${locus}/wtax4 > ./w_model_${locus}/wtax1
  perl ${SCRIPTS}/protaxscripts/thintaxonomy.pl 2 ./w_model_${locus}/wtax4 > ./w_model_${locus}/wtax2
  perl ${SCRIPTS}/protaxscripts/thintaxonomy.pl 3 ./w_model_${locus}/wtax4 > ./w_model_${locus}/wtax3
done
echo ""
echo "Step 2: generate training data"
for file in ${FILES[@]}
do
  locus=$(basename $file ".fa" | cut -f3 -d".")
  echo "Working on ${locus}"
  perl ${SCRIPTS}/protaxscripts/initialseqid2tax.pl ./w_model_${locus}/taxonomy ${taxon}.final_database.${locus}.fa > ./w_model_${locus}/seq2tax4
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL ${LEVEL}"
    perl ${SCRIPTS}/protaxscripts/makeseqid2tax.pl $LEVEL ./w_model_${locus}/seq2tax4 > ./w_model_${locus}/ref.wtax$LEVEL
    perl ${SCRIPTS}/protaxscripts/get1layer_reference_sequences_all.pl $LEVEL ./w_model_${locus}/wtax$LEVEL ./w_model_${locus}/ref.wtax$LEVEL ./w_model_${locus}/rseqs$LEVEL
    perl ${SCRIPTS}/protaxscripts/generate_training_data.pl ./w_model_${locus}/wtax$LEVEL ./w_model_${locus}/ref.wtax$LEVEL ./w_model_${locus}/rseqs$LEVEL 4500 1 no ./w_model_${locus}/train.w.level$LEVEL
    perl ${SCRIPTS}/protaxscripts/generate_unk_training_data.pl $LEVEL ./w_model_${locus}/wtax$LEVEL ./w_model_${locus}/ref.wtax$LEVEL ./w_model_${locus}/rseqs$LEVEL 500 1 no ./w_model_${locus}/train.w.unk$LEVEL
    cat ./w_model_${locus}/train.w.level$LEVEL ./w_model_${locus}/train.w.unk$LEVEL > ./w_model_${locus}/train.w$LEVEL
    cut -f6 -d" " ./w_model_${locus}/train.w$LEVEL | sort | uniq > ./w_model_${locus}/train.w${LEVEL}.id
  done
  cat ./w_model_${locus}/train.w[1,2,3,4].id | sort | uniq > ./w_model_${locus}/train.w.ids
  perl ${SCRIPTS}/protaxscripts/fastagrep.pl ./w_model_${locus}/train.w.ids ${taxon}.final_database.${locus}.fa > ./w_model_${locus}/training_${locus}.fa
done
echo ""
echo "Step 3: run LAST searches"
echo "This may take some time..."
for file in ${FILES[@]}
do
  locus=$(basename $file ".fa" | cut -f3 -d".")
  echo "Working on ${locus}"
  lastdb ./w_model_${locus}/lastref_${locus} ${taxon}.final_database.${locus}.fa
  lastal -T 1 -a 1 -f 0 -m 1000 ./w_model_${locus}/lastref_${locus} ./w_model_${locus}/training_${locus}.fa > ./w_model_${locus}/training.last
  perl ${SCRIPTS}/protaxscripts/last2sim.pl ./w_model_${locus}/training.last > ./w_model_${locus}/train.w.lastsim
  rm ./w_model_${locus}/training.last
done
echo ""
echo "Step 4: make x-matrices for input to R"
for file in ${FILES[@]}
do
  locus=$(basename $file ".fa" | cut -f3 -d".")
  echo "Working on ${locus}"
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL ${LEVEL}"
    perl ${SCRIPTS}/protaxscripts/create1layer_xdata4.pl ./w_model_${locus}/train.w$LEVEL ./w_model_${locus}/wtax$LEVEL ./w_model_${locus}/ref.wtax$LEVEL ./w_model_${locus}/rseqs$LEVEL ./w_model_${locus}/train.w.lastsim ./w_model_${locus}/train.w${LEVEL}.xdat 1
  done
done
echo ""
echo "Step 5: parameterise model"
echo "Running protax_training.R..."
Rscript --vanilla ${SCRIPTS}/weighted_protax_training.R ${SCRIPTS}/protaxscripts ${taxon}
echo ""

##### END
echo "End of train_weighted_protax.sh"
echo ""
end=`date +%s`
runtime=$((end-start))
echo "This took a total of `echo $runtime | awk '{printf "%.2f", $1/60}'` minutes (`echo $runtime | awk '{printf "%.2f", $1/3600}'` hours)."
echo ""
echo "Please select an mcmc iteration for each of the four levels for each marker (labelled ./w_model_${locus}/mcmc1a-d, ./w_model_${locus}/mcmc2a-d etc) based on the training plots (labelled ./w_model_${locus}/weighted_training_plot_${locus}_level1a_MCMC.pdf etc). Chains should be well-mixed and k-steps as close to 0.44 as possible. Relabel the selected model as ./w_model_${locus}/mcmc1 ./w_model_${locus}/mcmc2 etc."
echo ""
echo "Next step: Check model training with check_protax_training.sh"
echo "Then analyse real data. Either:"
echo "  - Process raw data with read_preprocessing.sh (experimental design must follow that described in the paper) and classify the output with weighted_protax_classify.sh"
echo "  - Classify OTUs with weighted_protax_classify_otus.sh"
echo ""
echo "Have a nice day :-)"
echo ""
