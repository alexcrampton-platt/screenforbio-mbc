#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a read_preprocessing.`date +%Y-%m-%d`.log)
exec 2> >(tee -a read_preprocessing.`date +%Y-%m-%d`.log >&2)

##### INFO

# read_preprocessing.sh

# for preprocessing PE illumina reads with twin tags (demultiplex, trim primers, QC, merge pairs)

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash read_preprocessing.sh pathToData MiSeq RunName readLength tagLength SampleSheet loci sampleTagDir
# where:
# pathToData is the path to the raw /Data/Intensities/BaseCalls directory
# MiSeq is the MiSeq maching name (e.g. M01108), the first 6 characters in the FASTQ read headers after '@'
# RunName should be an appropriate short label to be used as a prefix for output files and folder
# readLength is the read length
# tagLength is the tag length
# SampleSheet is a text file with information on the sample names and the sequencing (plate) tags - see bundled example for format.
# loci codes the combination of loci to search for as a binary string with each locus coded as either present (1) or absent (0) in the order: 12S,16S,CytB,COI. e.g. 16S only would be 0100, 12S and CytB would be 1010.
# sampleTagDir is a directory containing a text file with a list of sample names and the sample tags for each plate tag. file name format: PlateLabel.txt (PlateLabel should match Sample_ID in SampleSheet.csv extactly). file contents format: sample_name fwdTag revTag

##### WELCOME
if [ "$#" -ne 8 -o ! -r ${6} ]
then
	echo ""
	echo "You are trying to use read_preprocessing.sh but have not supplied enough/correct information. Please check the following:"
	echo ""
	echo "usage: bash read_preprocessing.sh pathToData MiSeq RunName readLength tagLength SampleSheet loci sampleTagDir"
	echo "where:"
	echo "pathToData is the path to the raw /Data/Intensities/BaseCalls directory"
	echo "MiSeq is the MiSeq maching name (e.g. M01108), the first 6 characters in the FASTQ read headers after '@'"
	echo "RunName should be an appropriate short label to be used as a prefix for output files and folder"
	echo "readLength is the read length"
	echo "tagLength is the tag length"
	echo "SampleSheet is a text file with information on the sample names and the sequencing (plate) tags - see bundled example for format."
	echo "loci codes the combination of loci to search for as a binary string with each locus coded as either present (1) or absent (0) in the order: 12S,16S,CytB,COI. e.g. 16S only would be 0100, 12S and CytB would be 1010"
	echo "sampleTagDir is a directory containing a text file with a list of sample names and the sample tags for each plate. file name format: PlateLabel.txt (PlateLabel should match Sample_ID in SampleSheet.csv extactly). file contents format: sample_name fwdTag revTag"
	echo ""
	exit 1
else
	echo ""
	echo "Welcome to read pre-processing for twin tagged Illumina amplicons"
	echo ""
	start=`date +%s`
	echo "The time now is $(date)."
	echo "Information from you:"
	echo "Raw data is in ${1}"
	echo "MiSeq name is ${2}"
	echo "File prefix and output directory will be ${3}"
	echo "Read length is ${4}, tag length is ${5}"
	echo "SampleSheet with plate tags is ${6}"
	echo "Primer sequences are in ${7}"
	echo "Sample tags are in ${8}"
fi

##### PARAMETERS
RAWDIR=${1}
MISEQ=${2}
PREFIX=${3}
READLEN=${4}
TAGLEN=${5}
PLATE=${6}
LOCI=${7}
SAMPLEDIR=${8}
DATE=`date +%Y-%m-%d`

SRNA=$(echo ${LOCI} | cut -c1)
LRNA=$(echo ${LOCI} | cut -c2)
CYTB=$(echo ${LOCI} | cut -c3)
COX1=$(echo ${LOCI} | cut -c4)

##### FUNCTIONS
function clip_12S {
	echo "	Clipping 12S primer sequences..."
	LOCUS=(12S)
	cutadapt -a AAAAAGCTTCAAACTGGGATTAGATACCCCACTAT...ACACACCGCCCGTCACCCTCTGCAGTCA$ --minimum-length 350 -o ./${PREFIX}/data/primerclip/${sample}.${LOCUS}.fq --discard-untrimmed ./${PREFIX}/data/merge/${sample}.merge.fq
	echo ""
	echo "	Done."
	echo ""
}
function clip_16S {
	echo "	Clipping 16S primer sequences..."
	LOCUS=(16S)
	cutadapt -a CGGTTGGGGTGACCTCGGA...AGTTACCCTAGGGATAACAGC$ --minimum-length 83 -o ./${PREFIX}/data/primerclip/${sample}.${LOCUS}.fq --discard-untrimmed ./${PREFIX}/data/merge/${sample}.merge.fq
	echo ""
	echo "	Done."
	echo ""
}
function clip_CytB {
	echo "	Clipping CytB primer sequences..."
	LOCUS=(CytB)
	cutadapt -a AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA...TGAGGACAAATATCATTCTGAGGGGCTGCAGTTT$ --minimum-length 271 -o ./${PREFIX}/data/primerclip/${sample}.${LOCUS}.fq --discard-untrimmed ./${PREFIX}/data/merge/${sample}.merge.fq
	echo ""
	echo "	Done."
	echo ""
}
function clip_COI {
	echo "	Clipping COI primer sequences..."
	LOCUS=(COI)
	cutadapt -a TCCACTAATCACAARGATATTGGTAC...ATAATCGGAGCCCCTGATA$ --minimum-length 200 -o ./${PREFIX}/data/primerclip/${sample}.${LOCUS}.fq --discard-untrimmed ./${PREFIX}/data/merge/${sample}.merge.fq
	echo ""
	echo "	Done."
	echo ""
}
function filter {
	usearch -fastq_filter ./${PREFIX}/data/primerclip/${sample}.${LOCUS}.fq -fastqout ./${PREFIX}/data/filter/${sample}.${LOCUS}.filter.fq -fastq_maxee 0.5 -threads 4
	echo ""
}
function derep {
	usearch -derep_fulllength ./${PREFIX}/data/filter/${sample}.${LOCUS}.filter.fq -fastqout ./${PREFIX}/data/derep/${sample}.${LOCUS}.filter.derep.fq -sizeout -strand both -minuniquesize 2 -relabel ${PREFIX}.${LOCUS}.${sample}_ -threads 4
	echo ""
}

##### MAIN
echo ""
echo "Start of read_preprocessing.sh"
echo ""
starting=`date +%s`
echo "The time now is $(date)"
echo ""
echo "Step 1: Setup directory"
echo "Generating output directory ./${PREFIX}..."
mkdir ./${PREFIX}
echo "Generating subdirectory structure..."
# make some subdirectories as needed
mkdir -p ./${PREFIX}/data/plates
mkdir -p ./${PREFIX}/data/samples
mkdir -p ./${PREFIX}/data/merge
mkdir -p ./${PREFIX}/data/primerclip
mkdir -p ./${PREFIX}/data/filter
mkdir -p ./${PREFIX}/data/derep
echo ""
echo "Output directory set up, now starting pre-processing."
echo ""
echo "Step 2: Base calls to FASTQ per plate"
echo "Running bcl2fastq, this may take some time..."
# using UMI setting
bcl2fastq --input-dir ${RAWDIR} --output-dir ./${PREFIX}/data/plates --barcode-mismatches 1 --with-failed-reads --minimum-trimmed-read-length ${READLEN} --sample-sheet ${PLATE} --loading-threads 4 --demultiplexing-threads 16 --writing-threads 4
echo ""
echo "Done. Compressed FASTQ files (R1 and R2 per plate) are in ./${PREFIX}/data/raw"
echo ""
echo "Checking file names..."
echo "Files are:"
file_path=($(find ./${PREFIX}/data/plates/${PREFIX} -mindepth 1 -maxdepth 1 -type f -name "*.f*q*"))
printf '%s\n' "${file_path[@]}"
echo ""
# Need to simplify names
echo "Renaming..."
for file in ./${PREFIX}/data/plates/${PREFIX}/*.f*q*
do mv "$file" "$(echo $file | sed 's/\(_S[0-9]\+_L001\)\(_R[12]\)\(_001\)/\2/g')"
done
echo "New names are:"
file_name=($(find ./${PREFIX}/data/plates/${PREFIX} -mindepth 1 -maxdepth 1 -type f -name "*.f*q*"))
printf '%s\n' "${file_name[@]}"
echo ""
echo "Plate names are:"
i=0
for name in ${file_name[@]}
do
	touch tmp
	plate_name[$i]="$(basename "$name" .fastq.gz)"
	printf '%s\n' "${plate_name[$i]}" | sed 's/_R[12]//g' - >> tmp
done
plate_name=($(cut -f 1 tmp | sort -u))
rm tmp
printf '%s\n' "${plate_name[@]}"
echo ""
echo "There are ${#plate_name[@]} plates to process."
echo ""
echo "Step 3: Assign to sample"
echo "Using AdapterRemoval to split reads from each plate into samples..."
for plate in ${plate_name[@]}
do
	echo "Processing $plate..."
	AdapterRemoval --file1 ./${PREFIX}/data/plates/${PREFIX}/${plate}_R1.fastq.gz --file2 ./${PREFIX}/data/plates/${PREFIX}/${plate}_R2.fastq.gz --basename ./${PREFIX}/data/samples/${plate} --barcode-list ${SAMPLEDIR}/${plate}.txt --barcode-mm-r1 1 --barcode-mm-r2 1 --maxns 0 --threads 4
done
echo ""
echo "Done."
echo ""
echo "Retrieving bad tag combinations in case these are useful for troubleshooting..."
for plate in ${plate_name[@]}
do
	grep "@"${MISEQ} ./${PREFIX}/data/samples/${plate}.unidentified_1 | cut -d ':' -f 8 | cut -d ' ' -f 1 | sort | uniq -c | sort -k1,1rn > ./${PREFIX}/data/samples/${plate}_bad_tag_combos.txt
	unID=$(grep -c "@"${MISEQ} ./${PREFIX}/data/samples/${plate}.unidentified_1)
	bad=$(grep -c "+" ./${PREFIX}/data/samples/${plate}_bad_tag_combos.txt)
	echo "${plate} has $bad bad tag combinations affecting $unID pairs."
done
echo ""
echo "Done. Bad combinations listed in ./${PREFIX}/data/samples/[plate]_bad_tag_combos.txt."
echo ""
# Rename files
echo "Renaming demultiplexed files to [sample].R[12].fq..."
for file in ./${PREFIX}/data/samples/*pair[12].truncated
do mv "$file" "$(echo $file | sed 's/\(.pair\)\([12]\)\(.truncated\)/.R\2.fq/g')"
done
echo ""
echo "Getting sample names..."
sample_path=($(find ./${PREFIX}/data/samples/ -mindepth 1 -maxdepth 1 -type f -name "*.R1.fq"))
for sample in ${sample_path[@]}
do
	touch tmp
	sample_name[$i]="$(basename "$sample" .R1.fq)"
	printf '%s\n' "${sample_name[$i]}" >> tmp
done
sample_name=($(cut -f 1 tmp))
rm tmp
printf '%s\n' "${sample_name[@]}"
echo ""
echo "Done. There are ${#sample_name[@]} samples to process further."
echo ""
echo "Step 4: Merge pairs"
echo "Merging pairs with usearch"
for sample in ${sample_name[@]}
do
	echo "Processing sample ${sample}..."
	if [ -s ./${PREFIX}/data/samples/${sample}.R1.fq ]
	then
		pairs=$(grep -c "@"${MISEQ} ./${PREFIX}/data/samples/${sample}.R1.fq)
		echo "	Sample has $pairs pairs"
		if [ $pairs -ge 1000 ]
		then
			usearch -fastq_mergepairs ./${PREFIX}/data/samples/${sample}.R1.fq -reverse ./${PREFIX}/data/samples/${sample}.R2.fq -fastqout ./${PREFIX}/data/merge/${sample}.merge.fq -fastq_minovlen 50 -fastq_maxdiffpct 20 -fastq_maxdiffs 20 -threads 4
			echo ""
		else
			echo "	Too few sequences (<1000) to process. Skipping..."
			echo ""
		fi
	else
		echo "	Demultiplexed reads for sample do not exist. Skipping..."
		echo ""
	fi
done
echo ""
echo "Done."
echo ""
echo "Step 5: Assign to locus and trim primers"
echo "Using cutadapt to determine locus and trim primers"
for sample in ${sample_name[@]}
do
	echo "Processing sample ${sample}..."
	if [ -f ./${PREFIX}/data/merge/${sample}.merge.fq ] && [ -s ./${PREFIX}/data/merge/${sample}.merge.fq ]
	then
		pairs=$(grep -c "@"${MISEQ} ./${PREFIX}/data/merge/${sample}.merge.fq)
		echo "	Sample has $pairs pairs"
		if [ $pairs -ge 500 ]
		then
			if [ ${SRNA} == 1 ]
			then
				echo "...12S present"
				clip_12S
			else
				echo "...12S absent"
			fi
			if [ ${LRNA} == 1 ]
			then
				echo "...16S present"
				clip_16S
			else
				echo "...16S absent"
			fi
			if [ ${CYTB} == 1 ]
			then
				echo "...CytB present"
				clip_CytB
			else
				echo "...CytB absent"
			fi
			if [ ${COX1} == 1 ]
			then
				echo "...COI present"
				clip_COI
			else
				echo "...COI absent"
			fi
		else
			echo "	Too few sequences (<500) to process. Skipping..."
			echo ""
		fi
	else
		echo "	No reads found. Skipping..."
	fi
done
echo ""
echo "Done."
echo ""
echo "Step 6: Quality filter"
echo "Filtering merged reads by expected error rate (must be ≤0.5) with usearch"
for sample in ${sample_name[@]}
do
	echo "Processing sample ${sample}..."
	if [ ${SRNA} == 1 ]
	then
		if [ -f ./${PREFIX}/data/primerclip/${sample}.12S.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.12S.fq ]
		then
			LOCUS=(12S)
			echo "	...filtering 12S"
			filter
		else
			echo "	...12S reads from previous step (primer clip) unavailable."
			echo ""
		fi
	else
		echo "...12S absent from run"
	fi
	if [ ${LRNA} == 1 ]
	then
		if [ -f ./${PREFIX}/data/primerclip/${sample}.16S.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.16S.fq ]
		then
			LOCUS=(16S)
			echo "	...filtering 16S"
			filter
		else
			echo "	...16S reads from previous step (primer clip) unavailable."
			echo ""
		fi
	else
		echo "...16S absent from run"
	fi
	if [ ${CYTB} == 1 ]
	then
		if [ -f ./${PREFIX}/data/primerclip/${sample}.CytB.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.CytB.fq ]
		then
			LOCUS=(CytB)
			echo "	...filtering CytB"
			filter
		else
			echo "	...CytB reads from previous step (primer clip) unavailable."
			echo ""
		fi
	else
		echo "...CytB absent from run"
	fi
	if [ ${COX1} == 1 ]
	then
		if [ -f ./${PREFIX}/data/primerclip/${sample}.COI.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.COI.fq ]
		then
			LOCUS=(COI)
			echo "	...filtering COI"
			filter
		else
			echo "	...COI reads from previous step (primer clip) unavailable."
			echo ""
		fi
	else
		echo "...COI absent from run"
	fi
done
echo ""
echo "Done."
echo ""
echo "Step 7: Dereplicate"
echo "Dereplicating filtered reads with usearch"
for sample in ${sample_name[@]}
do
	echo "Processing sample ${sample}..."
	if [ ${SRNA} == 1 ]
	then
		if [ -f ./${PREFIX}/data/filter/${sample}.12S.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.12S.filter.fq ]
		then
			LOCUS=(12S)
			echo "	...dereplicating 12S"
			derep
		else
			echo "	...12S reads from previous step (quality filter) unavailable."
			echo ""
		fi
	else
		echo "...12S absent from run"
	fi
	if [ ${LRNA} == 1 ]
	then
		if [ -f ./${PREFIX}/data/filter/${sample}.16S.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.16S.filter.fq ]
		then
			LOCUS=(16S)
			echo "	...dereplicating 16S"
			derep
		else
			echo "	...16S reads from previous step (quality filter) unavailable."
			echo ""
		fi
	else
		echo "...16S absent from run"
	fi
	if [ ${CYTB} == 1 ]
	then
		if [ -f ./${PREFIX}/data/filter/${sample}.CytB.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.CytB.filter.fq ]
		then
			LOCUS=(CytB)
			echo "	...dereplicating CytB"
			derep
		else
			echo "	...CytB reads from previous step (quality filter) unavailable."
			echo ""
		fi
	else
		echo "...CytB absent from run"
	fi
	if [ ${COX1} == 1 ]
	then
		if [ -f ./${PREFIX}/data/filter/${sample}.COI.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.COI.filter.fq ]
		then
			LOCUS=(COI)
			echo "	...dereplicating COI"
			derep
		else
			echo "	...COI reads from previous step (quality filter) unavailable."
			echo ""
		fi
	else
		echo "...COI absent from run"
	fi
done
echo ""
echo "Done."
echo ""
echo "Step 8: Print results table"
echo "Generate empty results table..."
touch ./${PREFIX}/pre-processing_results.${DATE}.txt
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Sample" "Raw_reads" "Merged_pairs" "12S_clipped" "16S_clipped" "CytB_clipped" "COI_clipped" "12S_filtered" "16S_filtered" "CytB_filtered" "COI_filtered" "12S_dereplicated" "16S_dereplicated" "CytB_dereplicateed" "COI_dereplicated" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
for sample in ${sample_name[@]}
do
	echo "Processing sample ${sample}..."
	echo "	...Raw reads..."
	if [ -f ./${PREFIX}/data/samples/${sample}.R1.fq ] && [ -s ./${PREFIX}/data/samples/${sample}.R1.fq ]
	then
		raw=$(grep -c "@"${MISEQ} ./${PREFIX}/data/samples/${sample}.R1.fq)
		printf "%s\t%s\t" "${sample}" "$raw" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t%s\t" "${sample}" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	echo "	...Merged pairs..."
	if [ -f ./${PREFIX}/data/merge/${sample}.merge.fq ] && [ -s ./${PREFIX}/data/merge/${sample}.merge.fq ]
	then
		merge=$(grep -c "@"${MISEQ} ./${PREFIX}/data/merge/${sample}.merge.fq)
		printf "%s\t" "$merge" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	echo "	...Primers clipped..."
	if [ -f ./${PREFIX}/data/primerclip/${sample}.12S.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.12S.fq ]
	then
		clip_srna=$(grep -c "@"${MISEQ} ./${PREFIX}/data/primerclip/${sample}.12S.fq)
		printf "%s\t" "$clip_srna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/primerclip/${sample}.16S.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.16S.fq ]
	then
		clip_lrna=$(grep -c "@"${MISEQ} ./${PREFIX}/data/primerclip/${sample}.16S.fq)
		printf "%s\t" "$clip_lrna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/primerclip/${sample}.CytB.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.CytB.fq ]
	then
		clip_cytb=$(grep -c "@"${MISEQ} ./${PREFIX}/data/primerclip/${sample}.CytB.fq)
		printf "%s\t" "$clip_cytb" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/primerclip/${sample}.COI.fq ] && [ -s ./${PREFIX}/data/primerclip/${sample}.COI.fq ]
	then
		clip_cox1=$(grep -c "@"${MISEQ} ./${PREFIX}/data/primerclip/${sample}.COI.fq)
		printf "%s\t" "$clip_cox1" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	echo "	...Filtered..."
	if [ -f ./${PREFIX}/data/filter/${sample}.12S.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.12S.filter.fq ]
	then
		filter_srna=$(grep -c "@"${MISEQ} ./${PREFIX}/data/filter/${sample}.12S.filter.fq)
		printf "%s\t" "$filter_srna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/filter/${sample}.16S.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.16S.filter.fq ]
	then
		filter_lrna=$(grep -c "@"${MISEQ} ./${PREFIX}/data/filter/${sample}.16S.filter.fq)
		printf "%s\t" "$filter_lrna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/filter/${sample}.CytB.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.CytB.filter.fq ]
	then
		filter_cytb=$(grep -c "@"${MISEQ} ./${PREFIX}/data/filter/${sample}.CytB.filter.fq)
		printf "%s\t" "$filter_cytb" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/filter/${sample}.COI.filter.fq ] && [ -s ./${PREFIX}/data/filter/${sample}.COI.filter.fq ]
	then
		filter_cox1=$(grep -c "@"${MISEQ} ./${PREFIX}/data/filter/${sample}.COI.filter.fq)
		printf "%s\t" "$filter_cox1" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	echo "	...Dereplicated..."
	if [ -f ./${PREFIX}/data/derep/${sample}.12S.filter.derep.fq ] && [ -s ./${PREFIX}/data/derep/${sample}.12S.filter.derep.fq ]
	then
		derep_srna=$(grep -c "^@p" ./${PREFIX}/data/derep/${sample}.12S.filter.derep.fq)
		printf "%s\t" "$derep_srna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/derep/${sample}.16S.filter.derep.fq ] && [ -s ./${PREFIX}/data/derep/${sample}.16S.filter.derep.fq ]
	then
		derep_lrna=$(grep -c "^@p" ./${PREFIX}/data/derep/${sample}.16S.filter.derep.fq)
		printf "%s\t" "$derep_lrna" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/derep/${sample}.CytB.filter.derep.fq ] && [ -s ./${PREFIX}/data/derep/${sample}.CytB.filter.derep.fq ]
	then
		derep_cytb=$(grep -c "^@p" ./${PREFIX}/data/derep/${sample}.CytB.filter.derep.fq)
		printf "%s\t" "$derep_cytb" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\t" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	if [ -f ./${PREFIX}/data/derep/${sample}.COI.filter.derep.fq ] && [ -s ./${PREFIX}/data/derep/${sample}.COI.filter.derep.fq ]
	then
		derep_cox1=$(grep -c "^@p" ./${PREFIX}/data/derep/${sample}.COI.filter.derep.fq)
		printf "%s\n" "$derep_cox1" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	else
		printf "%s\n" "0" >> ./${PREFIX}/pre-processing_results.${DATE}.txt
	fi
	echo "	...Done."
	echo ""
done
echo "Done."
echo ""

##### END
echo "End of read_processing.sh"
echo ""
ending=`date +%s`
echo "The time now is $(date)"
runningtime=$((ending-starting))
echo ""
echo "Read processing took `echo $runningtime | awk '{printf "%.2f", $1/60}'` minutes (`echo $runningtime | awk '{printf "%.2f", $1/3600}'` hours) to complete."
echo ""
echo "Next step: Taxonomic assignment with PROTAX. Choice of analysis depends on the models to be used:"
echo "	Unweighted - protax_classify.sh"
echo "	Weighted - weighted_protax_classify.sh"
echo ""
echo "Have a nice day :-)"
echo ""
