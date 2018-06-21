#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a protax_classify.`date +%Y-%m-%d`.log)
exec 2> >(tee -a protax_classify.`date +%Y-%m-%d`.log >&2)

##### INFO

# protax_classify.sh

# for classifying batches of pre-preprocessed reads using a pre-trained unweighted protax model and cleaned reference database

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash protax_classify.sh pathToData locus protaxdir screenforbio outdir
# where:
# pathToData is the path to a directory containing FASTQ files to be processed one by one
# locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus per run, run script once for each.
# protaxdir is the path to a directory containing protax models and clean databases for all 4 loci
# screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
# outdir is the name to give an output directory (inside current)

##### WELCOME
if [ "$#" -ne 5 ]
then
	echo ""
	echo "You are trying to use protax_classify.sh but have not provided enough information."
	echo "Please check the following:"
	echo ""
	echo "usage: bash protax_classify.sh pathToData locus protaxdir outdir"
	echo "where:"
	echo "pathToData is the path to a directory containing FASTQ files to be processed one by one"
  echo "locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus per run, run script once for each."
  echo "protaxdir is the path to a directory containing protax models and clean databases for all 4 loci (one subdirectory per locus) plus required scripts"
	echo "screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)"
	echo "outdir is the name to give an output directory (inside current)"
	echo ""
	exit 1
else
	echo "Welcome to protax_classify.sh: batch classification of FASTQ files with an unweighted PROTAX model"
	echo ""
	echo "The time now is $(date)."
	echo "Information from you:"
	echo "FASTQ files to be processed are in ${1}"
  echo "Protax models and reference databases are in ${3}"
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
echo "Step 1: Get names of input files"
echo ""
echo "File names are:"
file_name=($(find ${DATA} -mindepth 1 -maxdepth 1 -type f -name "*.${LOCUS}.filter.derep.fq"))
printf '%s\n' "${file_name[@]}"
echo ""
echo "There are ${#file_name[@]} files to process."
echo ""

echo "Step 2: Make output directory"
echo ""
mkdir ./${OUTDIR}_${LOCUS}
echo "Output will be in: ./${OUTDIR}_${LOCUS}"
echo ""

echo "Step 3: Run protax on each input file"
echo ""
echo "This may take a while..."
echo ""
start=`date +%s`
for file in ${file_name[@]}
do
	if [ -s ${file} ]
	then
		echo "Processing file: ${file}"
		starting=`date +%s`
		echo "	Starting at $(date)."
		#get sample name
		sample_name="$(basename "$file" .filter.derep.fq)"
	  #make output directory for sample
	  mkdir ./${OUTDIR}_${LOCUS}/${sample_name}/
	  #convert fastq to fasta
		echo "	Converting FASTQ to FASTA..."
	  usearch -fastq_filter ${file} -fastaout ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.fa
	  #run LAST search
		echo "	Running LAST search..."
	  lastal -T 1 -a 1 -f 0 -m 1000 ${MODELS}/model_${LOCUS}/lastref_${LOCUS} ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.fa > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.last
	  #convert LAST result to similarity
		echo "	Convert LAST result to sequence similarity..."
	  perl ${SCRIPTS}/protaxscripts/last2sim.pl ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.last > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.lastsim
	  #clean up
	  rm ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.last
	  #get IDs
		echo "	Get sequence IDs..."
	  grep '^>' ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.fa | cut -c2- > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.ids
	  #generate base log probability file
		echo "	Generate base log probability file..."
	  perl ${SCRIPTS}/protaxscripts/testsample2init.pl ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.ids > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.0.logprob
	  #classify sequences at each taxonomic level (order, family, genus, species)
		echo "	Classifying at..."
	  for LEVEL in 1 2 3 4
	  do
	    echo "		LEVEL $LEVEL"
	    PREVLEVEL=$((LEVEL-1))
	    IFILE=./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.${PREVLEVEL}.logprob
	    OFILE=./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.${LEVEL}.logprob
	    perl ${SCRIPTS}/protaxscripts/classify4.pl $IFILE ${MODELS}/model_${LOCUS}/tax$LEVEL ${MODELS}/model_${LOCUS}/ref.tax$LEVEL ${MODELS}/model_${LOCUS}/rseqs$LEVEL ${MODELS}/model_${LOCUS}/mcmc$LEVEL map ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.lastsim 0 .05 $OFILE 1
	  done
	  #add taxonomic info
		echo "	Adding taxonomy..."
	  perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.0.logprob > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.class_probs
	  perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.1.logprob > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.order_probs
	  perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.2.logprob > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.family_probs
	  perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.3.logprob > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.genus_probs
	  perl ${SCRIPTS}/protaxscripts/add_taxonomy_info.pl ${MODELS}/model_${LOCUS}/taxonomy ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.4.logprob > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.species_probs
		#add sequence similarity of assigned species/genus (for unassigned reads takes best matching sequence)
		for id in $(cut -f1 ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.species_probs | sort -u)
		do
			sp=$(grep ${id} ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.species_probs | cut -f5 | awk 'BEGIN{FS=","}{print $3 "_" $4}')
			if grep -q "_unk$" <(echo ${sp})
			then
				gen=$(sed 's/_unk//' <(echo ${sp}))
				#find gen
				grep ${id} ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.lastsim | grep ${gen} | awk -v max=0 '{if($3>max){want=$0;max=$3}}END{print want}' >> ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.bestsim
			else
				#find sp (for unassigned reads will pick up best matching sequence)
				grep ${id} ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.lastsim | grep ${sp} | awk -v max=0 '{if($3>max){want=$0;max=$3}}END{print want}' >> ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.bestsim
			fi
		done
		join -j 1 -a 1 -o 1.1,1.2,1.3,1.4,1.5,2.3,2.2 <(sort -k1,1 ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.species_probs) <(sort -k1,1 ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.bestsim) > ./${OUTDIR}_${LOCUS}/${sample_name}/${sample_name}.species_probs_sim
		#finish iteration with message
		echo ""
		echo "	Finished at $(date)."
		ending=`date +%s`
		runningtime=$((ending-starting))
		echo "	Processing of ${file} took `echo $runningtime | awk '{printf "%.2f", $1/60}'` minutes"
		echo ""
	else
		echo "${file} empty, skipping..."
		echo ""
	fi
done
echo ""

##### END
echo "protax_classify complete for all input files"
end=`date +%s`
runtime=$((end-start))
echo "This took a total of `echo $runtime | awk '{printf "%.2f", $1/60}'` minutes (`echo $runtime | awk '{printf "%.2f", $1/3600}'` hours)."
echo ""
echo "Result files are in ./${OUTDIR}_${LOCUS}/sample_name"
echo "Classification for each read at each taxonomic level (species, genus, family, order) in files sample_name.level_probs"
echo "Headers are:"
echo "	queryID	taxID	log(probability)	level	taxon"
echo "Additionally, the best matching hit (for assigned species/genus where available) found with LAST is appended to sample_name.species_probs in sample_name.species_probs_sim"
echo "Headers are:"
echo "	queryID	taxID	log(probability)	level	taxon	bestHit_similarity	bestHit"
echo ""
echo "Have a nice day :-)"
echo ""
