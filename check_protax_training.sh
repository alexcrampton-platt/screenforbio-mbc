#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a check_protax_training.`date +%Y-%m-%d`.log)
exec 2> >(tee -a check_protax_training.`date +%Y-%m-%d`.log >&2)

##### INFO

# check_protax_training.sh

# for generating classification accuracy plots for training data with selected PROTAX model

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash check_protax_training.sh modeldir taxon locus screenforbio
# where:
# modeldir is the path to a directory containing the protax model to be checked
# taxon is the taxon for which the model was generated (used for labelling only)
# locus is the locus for which the model was generated (used for labelling only)
# screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)

##### WELCOME
if [ "$#" -ne 4 ]
then
	echo ""
	echo "You are trying to use check_protax_training.sh but have not provided enough information."
	echo "Please check the following:"
	echo ""
	echo "Usage: bash check_protax_training.sh modeldir taxon locus screenforbio"
	echo "Where:"
	echo "modeldir is the path to a directory containing the protax model to be checked"
	echo "taxon is the taxon for which the model was generated (used for labelling only)"
  echo "locus is the locus for which the model was generated (used for labelling only)"
	echo "screenforbio is the path to the screenforbio-mbc folder (must contain subdirectory protaxscripts)"
	echo ""
	exit 1
else
	echo ""
	echo "Welcome to check_protax_training.sh: generating classification accuracy plots for specified PROTAX model"
	echo ""
	echo "The time now is $(date)."
	echo ""
	echo "Generating plots for model in ${1}"
	echo ""
fi

##### PARAMETERS
MODEL=${1}
TAXON=${2}
LOCUS=${3}
SCRIPTS=${4}

##### MAIN
echo "Check if model is weighted or unweighted"
echo ""
start=`date +%s`
if [ $(basename ${MODEL} | cut -f1 -d"_") = "w" ]
then
  modtype=$(echo "weighted")
else
  modtype=$(echo "unweighted")
fi
echo "Model is ${modtype}"
echo ""
echo "Classifying training data with model"
echo ""
mkdir ${MODEL}/checktrain
if [ ${modtype} = "unweighted" ]
then
  perl ${SCRIPTS}/protaxscripts/trainsample2init.pl ${MODEL}/train4 > ${MODEL}/checktrain/query0.logprob
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL $LEVEL"
    PREVLEVEL=$((LEVEL-1))
    IFILE=${MODEL}/checktrain/query${PREVLEVEL}.logprob
    OFILE=${MODEL}/checktrain/query${LEVEL}.logprob
    perl ${SCRIPTS}/protaxscripts/trainclassify4.pl $IFILE ${MODEL}/tax$LEVEL ${MODEL}/ref.tax$LEVEL ${MODEL}/rseqs$LEVEL ${MODEL}/mcmc$LEVEL map ${MODEL}/train.lastsim 0 .01 $OFILE 1
  done
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL $LEVEL"
    perl ${SCRIPTS}/protaxscripts/trainsample2correct.pl $LEVEL ${MODEL}/tax4 ${MODEL}/train4 > ${MODEL}/checktrain/query${LEVEL}.tax
    perl ${SCRIPTS}/protaxscripts/trainsample2addcor.pl ${MODEL}/checktrain/query${LEVEL}.logprob ${MODEL}/checktrain/query${LEVEL}.tax ${MODEL}/tax$LEVEL > ${MODEL}/checktrain/query${LEVEL}.cor
  done
else
  perl ${SCRIPTS}/protaxscripts/trainsample2init.pl ${MODEL}/train.w4 > ${MODEL}/checktrain/query0.logprob
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL $LEVEL"
    PREVLEVEL=$((LEVEL-1))
    IFILE=${MODEL}/checktrain/query${PREVLEVEL}.logprob
    OFILE=${MODEL}/checktrain/query${LEVEL}.logprob
    perl ${SCRIPTS}/protaxscripts/trainclassify4.pl $IFILE ${MODEL}/wtax$LEVEL ${MODEL}/ref.wtax$LEVEL ${MODEL}/rseqs$LEVEL ${MODEL}/w_mcmc$LEVEL map ${MODEL}/train.w.lastsim 0 .01 $OFILE 1
  done
  for LEVEL in 1 2 3 4
  do
    echo "LEVEL $LEVEL"
    perl ${SCRIPTS}/protaxscripts/trainsample2correct.pl $LEVEL ${MODEL}/wtax4 ${MODEL}/train.w4 > ${MODEL}/checktrain/query${LEVEL}.tax
    perl ${SCRIPTS}/protaxscripts/trainsample2addcor.pl ${MODEL}/checktrain/query${LEVEL}.logprob ${MODEL}/checktrain/query${LEVEL}.tax ${MODEL}/wtax$LEVEL > ${MODEL}/checktrain/query${LEVEL}.cor
  done
fi
echo ""
echo "Generating plots with training_plots.R..."
echo ""
Rscript --vanilla ${SCRIPTS}/training_plots.R ${MODEL}/checktrain ${TAXON} ${LOCUS} ${modtype}
echo ""
echo "...training_plots.R finished"
end=`date +%s`
runtime=$((end-start))
echo ""

##### END
echo "Model check took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
echo ""
echo "Plots can be found in ${MODEL}/checktrain/${modtype}_${TAXON}_${LOCUS}_biasaccuracy.pdf"
echo ""
echo "Next step: analyse real data. Either:"
echo "  - Process raw data with read_preprocessing.sh (experimental design must follow that described in the manuscript) and classify the output with protax_classify.sh or weighted_protax_classify.sh as appropriate"
echo "  - Classify OTUs with protax_classify_otus.sh or weighted_protax_classify_otus.sh as appropriate"
echo ""
echo "Have a nice day :-)"
echo ""
