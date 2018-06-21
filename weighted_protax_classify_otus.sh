#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a weighted_protax_classify_otus.`date +%Y-%m-%d`.log)
exec 2> >(tee -a weighted_protax_classify_otus.`date +%Y-%m-%d`.log >&2)

##### INFO

# weighted_protax_classify_otus.sh

# for classifying OTUs using a pre-trained weighted protax model and cleaned reference database

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash weighted_protax_classify_otus.sh otus locus protaxdir screenforbio outdir
# where:
# otus is the (path to) the OTU fasta to be processed (suffix should be ".fa")
# locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus to analyse, run script once for each.
# protaxdir is the path to a directory containing weighted protax models and clean databases for all 4 loci
# screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
# outdir is the name to give the output directory (inside current)

##### WELCOME
if [ "$#" -ne 5 ]
then
	echo ""
	echo "You are trying to use weighted_protax_classify_otus.sh but have not provided enough information."
	echo "Please check the following:"
	echo ""
	echo "usage: bash weighted_protax_classify_otus.sh otus locus protaxdir screenforbio outdir"
	echo "where:"
	echo "otus is the OTU FASTA to be processed (suffix should be ".fa")."
  echo "locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus to analyse, run script once for each."
  echo "protaxdir is the path to a directory containing weighted protax models and clean databases for all 4 loci (one subdirectory per locus)"
	echo "screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)"
	echo "outdir is the name to give the output directory (inside current)"
	echo ""
	exit 1
else
	echo "Welcome to weighted_protax_classify_otus.sh: classification of OTUs with a weighted PROTAX model"
	echo ""
	echo "The time now is $(date)."
	echo "Information from you:"
	echo "FASTA file of OTUs is ${1}"
	echo "Weighted protax models and reference databases are in ${3}"
	echo "PROTAX Perl scripts are in ${4}/protaxscripts"
	echo "Output directory will be ${5}"
	echo ""
fi

##### PARAMETERS
DATA=${1}
LOCUS=${2}
MODELS=${3}
SCRIPTS=${4}
OUTDIR=${5}
DATE=`date +%Y-%m-%d`

##### MAIN
echo "Step 1: Get names of input file"
echo ""
echo "File names is:"
file_name=($(basename ${DATA} ".fa"))
printf '%s\n' "${file_name[@]}"
echo ""

echo "Step 2: Make output directory"
echo ""
mkdir ./${OUTDIR}_${LOCUS}
echo "Output will be in: ./${OUTDIR}_${LOCUS}"
echo ""

echo "Step 3: Run protax"
echo ""
start=`date +%s`
#run LAST search
echo "	Running LAST search..."
lastal -T 1 -a 1 -f 0 -m 1000 ${MODELS}/w_model_${LOCUS}/lastref_${LOCUS} ${DATA} > ./${OUTDIR}_${LOCUS}/${file_name}.last
#convert LAST result to similarity
echo "	Convert LAST result to sequence similarity..."
perl ${SCRIPTS}/protaxscripts/last2sim.pl ./${OUTDIR}_${LOCUS}/${file_name}.last > ./${OUTDIR}_${LOCUS}/${file_name}.lastsim
#clean up
rm ./${OUTDIR}_${LOCUS}/${file_name}.last
#get IDs
echo "	Get sequence IDs..."
grep '^>' ${DATA} | cut -c2- > ./${OUTDIR}_${LOCUS}/${file_name}.ids
#generate base log probability file
echo "	Generate base log probability file..."
perl ${SCRIPTS}/protaxscripts/testsample2init.pl ./${OUTDIR}_${LOCUS}/${file_name}.ids > ./${OUTDIR}_${LOCUS}/${file_name}.0.w_logprob
#classify sequences at each taxonomic level (order, family, genus, species)
echo "	Classifying at..."
for LEVEL in 1 2 3 4
do
  echo "		LEVEL $LEVEL"
  PREVLEVEL=$((LEVEL-1))
  IFILE=./${OUTDIR}_${LOCUS}/${file_name}.${PREVLEVEL}.w_logprob
  OFILE=./${OUTDIR}_${LOCUS}/${file_name}.${LEVEL}.w_logprob
  perl ${SCRIPTS}/protaxscripts/classify4.pl $IFILE ${MODELS}/w_model_${LOCUS}/wtax$LEVEL ${MODELS}/w_model_${LOCUS}/ref.wtax$LEVEL ${MODELS}/w_model_${LOCUS}/rseqs$LEVEL ${MODELS}/w_model_${LOCUS}/w_mcmc$LEVEL map ./${OUTDIR}_${LOCUS}/${file_name}.lastsim 0 .05 $OFILE 1
done
#add taxonomic info
echo "	Adding taxonomy..."
perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/w_model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${file_name}.0.w_logprob > ./${OUTDIR}_${LOCUS}/${file_name}.w_class_probs
perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/w_model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${file_name}.1.w_logprob > ./${OUTDIR}_${LOCUS}/${file_name}.w_order_probs
perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/w_model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${file_name}.2.w_logprob > ./${OUTDIR}_${LOCUS}/${file_name}.w_family_probs
perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/w_model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${file_name}.3.w_logprob > ./${OUTDIR}_${LOCUS}/${file_name}.w_genus_probs
perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/w_model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${file_name}.4.w_logprob > ./${OUTDIR}_${LOCUS}/${file_name}.w_species_probs
#add sequence similarity of assigned species/genus (for unassigned OTUs takes best matching sequence)
for id in $(cut -f1 ./${OUTDIR}_${LOCUS}/${file_name}.w_species_probs | sort -u)
do
	sp=$(grep -w ${id} ./${OUTDIR}_${LOCUS}/${file_name}.w_species_probs | cut -f5 | awk 'BEGIN{FS=","}{print $3 "_" $4}')
	if grep -q "_unk$" <(echo ${sp})
	then
		gen=$(sed 's/_unk//' <(echo ${sp}))
		#find gen
		if grep -q -w ${id} ./${OUTDIR}_${LOCUS}/${file_name}.lastsim
		then
			grep -w ${id} ./${OUTDIR}_${LOCUS}/${file_name}.lastsim | grep ${gen} | awk -v max=0 '{if($3>max){want=$0;max=$3}}END{print want}' >> ./${OUTDIR}_${LOCUS}/${file_name}.bestsim
		else
			continue
		fi
	else
		#find sp (for unassigned OTUs will pick up best matching sequence)
		if grep -q -w ${id} ./${OUTDIR}_${LOCUS}/${file_name}.lastsim
		then
			grep -w ${id} ./${OUTDIR}_${LOCUS}/${file_name}.lastsim | grep ${sp} | awk -v max=0 '{if($3>max){want=$0;max=$3}}END{print want}' >> ./${OUTDIR}_${LOCUS}/${file_name}.bestsim
		else
			continue
		fi
	fi
done
join -j 1 -a 1 -o 1.1,1.2,1.3,1.4,1.5,2.3,2.2 <(sort -k1,1 ./${OUTDIR}_${LOCUS}/${file_name}.w_species_probs) <(sort -k1,1 ./${OUTDIR}_${LOCUS}/${file_name}.bestsim) > ./${OUTDIR}_${LOCUS}/${file_name}.w_species_probs_sim
echo ""

##### END
echo "weighted_protax_classify_otus.sh complete."
end=`date +%s`
runtime=$((end-start))
echo "This took a total of `echo $runtime | awk '{printf "%.2f", $1/60}'` minutes."
echo ""
echo "Results are in ./${OUTDIR}_${LOCUS}"
echo "Classification for each OTU at each taxonomic level (species, genus, family, order) in files ${file_name}.w_level_probs"
echo "Headers are:"
echo "	queryID	taxID	log(probability)	level	taxon"
echo "Additionally, the best matching hit (for assigned species/genus where available) found with LAST is appended to ${file_name}.w_species_probs in ${file_name}.w_species_probs_sim"
echo "Headers are:"
echo "	queryID	taxID	log(probability)	level	taxon	bestHit_similarity	bestHit"
echo ""
echo "Have a nice day :-)"
echo ""
