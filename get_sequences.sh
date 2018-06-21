#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a get_sequences_module-${3}.`date +%Y-%m-%d`.log)
exec 2> >(tee -a get_sequences_module-${3}.`date +%Y-%m-%d`.log >&2)

##### INFO

# get_sequences.sh

# for generating reference databases based primarily on MIDORI UNIQUE sequences with RDP headers (available from http://www.reference-midori.info/download.php#)

# written by Alex Crampton-Platt for Andreas Wilting (IZW, ScreenForBio project)

# usage: bash get_sequences.sh extras gap-fill module taxon screenforbio
# where:
# extras is 'yes' or 'no', indicating whether to add local FASTA format sequences. if 'yes', files must be in present directory labelled "extra_12S.fa", "extra_16S.fa", "extra_Cytb.fa", "extra_COI.fa", with headers in format Genus_species_uniqueID.
# gap-fill is 'no' or a tab-delimited text file of species names to be targeted for gap-filling from NCBI, in format Genus_species.
# module is 'one', 'two', 'three' or 'four' indicating whether the script is starting from scratch ('one'), restarting after checking the output of the mafft alignment ('two'), restarting after manual correction of failed taxonomy lookups ('three'), or restartinh after manual checks of SATIVA output ('four'). see end of module messages for any requirements for the next module."
# taxon is the taxon for which the taxonomy was downloaded with get_taxonomy.sh, e.g. Mammalia or Tetrapoda (all outputs should be in present directory).
# screenforbio is the path to the screenforbio-mbc directory

# assumes that the file are in current directory and names are in format MIDORI_UNIQUE_x.x_locus_RDP.fasta where x.x is version number.
# assumes collapsetypes_v4.6.pl is in screenforbio directory - available from https://sourceforge.net/projects/collapsetypes/

##### WELCOME
if [ "$#" -ne 5 ]
then
  echo ""
  echo "You are trying to use get_sequences.sh but have not supplied enough. Please check the following:"
  echo ""
  echo "usage: bash get_sequences.sh extras gap-fill module taxon screenforbio"
  echo "where:"
  echo "extras is 'yes' or 'no', indicating whether to add local FASTA format sequences. if 'yes', files must be in present directory labelled "extra_12S.fa", "extra_16S.fa", "extra_Cytb.fa", "extra_COI.fa", with headers in format Genus_species_uniqueID."
  echo "gap-fill is 'no' or a tab-delimited text file of species names to be targeted for gap-filling from NCBI, in format Genus_species."
  echo "module is 'one', 'two', 'three' or 'four' indicating whether the script is starting from scratch ('one'), restarting after checking the output of the mafft alignment ('two'), restarting after manual correction of failed taxonomy lookups ('three'), or restartinh after manual checks of SATIVA output ('four'). see end of module messages for any requirements for the next module."
  echo "taxon is the taxon for which the taxonomy was downloaded with get_taxonomy.sh, e.g. Mammalia or Tetrapoda (all outputs should be in present directory)."
  echo "screenforbio is the path to the screenforbio-mbc directory"
  echo ""
  exit 1
else
  echo ""
  echo "Welcome to get_sequences.sh: generating a curated reference sequence database"
  echo ""
fi

##### PARAMETERS
EXTRAS=${1}
GAPFILL=${2}
MODULE=${3}
TAXON=${4}
SCRIPTS=${5}
DATE=`date +%Y-%m-%d`

##### FUNCTIONS
function module_one {
  echo "Starting from scratch with raw FASTA files..."
  echo ""
  echo "  Looking for sequence FASTAs..."
  if ls MIDORI_UNIQUE_*.fasta 1> /dev/null 2>&1
  then
  FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_UNIQUE_*.fasta" | sed 's/\.\///g'))
  echo "  Done."
  else
    echo "  No matching FASTAs found."
    exit 1
  fi
  echo ""
  echo "Processing from raw FASTA to aligned amplicons..."
  echo ""
  for file in ${FILES[@]}
  do
    echo "Working on ${file}..."
    label=($(echo ${file} | sed 's/[0-9]\.[0-9]_//g' | sed 's/_RDP.fasta//g' | sed 's/UNIQUE_//g'))
    locus=($(echo ${label} | sed 's/MIDORI_//g'))
    #get ${TAXON} subset
    if [ ${TAXON} = "Tetrapoda" ]
    then
      grep -E "Amphibia|Aves|Crocodylia|Mammalia|Sphenodontia|Squamata|Testudines" ${file} | sed 's/>//g' | sed 's/ /_/g' > ${label}.headers.txt
    else
      if grep -q -m1 -w ${TAXON} ${file}
      then
        grep ${TAXON} ${file} | sed 's/>//g' | sed 's/ /_/g' > ${label}.headers.txt
      else
        echo "  No ${TAXON} sequences found in ${file}, skipping..."
        echo ""
        continue
      fi
    fi
    #strip headers to species only (part after last ;)
    sed 's/root;Eukaryota;[A-Za-z].\+;//g' ${label}.headers.txt > ${label}.species.txt
    #find non-binomial labels
    if grep -q -m1 -E "/|_x_|_[A-Z]|_[0-9]|-[0-9]|_\(|\." ${label}.species.txt
    then
      grep -E "/|_x_|_[A-Z]|_[0-9]|-[0-9]|_\(|\." ${label}.species.txt > ${label}.drop.txt
    else
      touch ${label}.drop.txt
    fi
    #get list of properly labelled sequences, reduce any subspecies to species
    if [ -s ${label}.drop.txt ]
    then
      tabtk isct -c ${label}.drop.txt ${label}.species.txt | sed 's/\(\t[A-Z][a-z].\+_[a-z].\+\)\(_[a-z].\+\)/\1/g' | sort -k1,1 -u | awk '{print $1 "\t" $2 "\t" $2 "_" $1}' > ${label}.keep.txt
    else
      awk '{print $1 "\t" $2 "\t" $2 "_" $1}' ${label}.species.txt > ${label}.keep.txt
    fi
    #pickup seqs
    seqtk subseq ${file} <(cut -f1 ${label}.keep.txt) | awk '{print $1}' > ${label}.keep.fa
    #remove duplicated seqs (same accession AND length - can have same accession w/different seqs in MIDORI...)
    awk '/^>/{f=!d[$1];d[$1]=1}f' <(seqkit replace -p '(.+)$' -r '{kv}' -k <(seqtk comp ${label}.keep.fa | awk '{print $1 "\t" $1 "_" $2}') ${label}.keep.fa --keep-key) > ${label}.keep_dedup.fa
    #rename
    seqkit replace -p '(.+)$' -r '{kv}' -k <(awk '{print $1 "\t" $3}' ${label}.keep.txt) <(sed 's/_.\+//g' ${label}.keep_dedup.fa) --keep-key --quiet > ${label}.keep_relabel.fa
    #cleanup
    rm ${label}.headers.txt
    rm ${label}.species.txt
    rm ${label}.drop.txt
    echo "  Finding amplicons..."
    #get amplicons
    if [ ${locus} = "srRNA" ]
    then
      echo "  Amplicon should be ${locus}..."
      # add extras if requested
      if [ ${EXTRAS} = "yes" ] && [ -s extra_12S.fa ]
      then
        echo "  Adding extra sequences..."
        # add extra sequences (make sure on single line and uppercase)
        cat ${label}.keep_relabel.fa <(awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' extra_12S.fa) | seqtk seq -l0 - > ${label}.raw.fa
      else
        echo "  No extras to add. Skipping..."
        cp ${label}.keep_relabel.fa ${label}.raw.fa
      fi
      if [ -s ${GAPFILL} ]
      then
        echo "  Checking for target species..."
        splist=($(cut -f1 ${GAPFILL}))
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"s-rRNA\" OR \"12S rRNA\" OR \"12S ribosomal RNA\"\)" | efilter -location mitochondrion -source genbank | efetch -format fasta >> target_species.12S.genbank.${DATE}.fa; done`
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"s-rRNA\" OR \"12S rRNA\" OR \"12S ribosomal RNA\"\)" | efilter -location mitochondrion -source refseq | efetch -format fasta >> target_species.12S.refseq.${DATE}.fa; done`
        cat target_species.12S.genbank.${DATE}.fa target_species.12S.refseq.${DATE}.fa > target_species.12S.ncbi.${DATE}.fa
        # edit NC_ to NC (and AC_ to AC)
        sed -i 's/\(>[A-Z][A-Z]\)\(_\)/\1/g' target_species.12S.ncbi.${DATE}.fa
        #discard dupicates based on ID
        awk '/^>/{f=!d[$1];d[$1]=1}f' target_species.12S.ncbi.${DATE}.fa > tmp
        #drop UNVERIFIED and any with {}s in header
        seqtk subseq tmp <(grep ">" tmp | sed 's/>//g' | sed '/UNVERIFIED/d' | sed '/{/d' | sed '/}/d') > target_species.12S.ncbi.${DATE}.fa
        rm tmp
        # get current headers and make a renaming file to convert to same format as protax database (Genus_species_accession)
        grep ">" target_species.12S.ncbi.${DATE}.fa | sed 's/>//g' > target_species.12S.ncbi.${DATE}.headers.txt
        cut -f1,2,3 -d" " target_species.12S.ncbi.${DATE}.headers.txt | sed 's/\.[0-9]/ /g' | awk '{print $1 "\t" $2 "_" $3 "_" $1}' > target_species.12S.ncbi.${DATE}.rename.txt
        # simplify current headers
        awk '{print $1}' target_species.12S.ncbi.${DATE}.fa | sed 's/\(^>[A-Z0-9].\+\)\(\.[0-9]\)/\1/g' > target_species.12S.ncbi.${DATE}.ids.fa
        # get list of accessions and check against list of previously seen accessions
        grep ">" target_species.12S.ncbi.${DATE}.ids.fa | sed 's/>//g' | sort -u | cut -f1 -d" " | cut -f1 -d"." > target_species.12S.ncbi.${DATE}.ids
        ids=($(cut -f1 target_species.12S.ncbi.${DATE}.ids))
        for id in ${ids[@]}
        do
          if grep -q ${id} ${label}.raw.fa
          then
            printf '%s\t%s\n' ${id} "true" >> target_species.12S.ncbi.${DATE}.ids_status.txt
          else
            printf '%s\t%s\n' ${id} "false" >> target_species.12S.ncbi.${DATE}.ids_status.txt
          fi
        done
        # select previously unseen sequences
        grep "false" target_species.12S.ncbi.${DATE}.ids_status.txt | cut -f1 > target_species.12S.ncbi.${DATE}.ids_novel.txt
        seqtk subseq target_species.12S.ncbi.${DATE}.ids.fa target_species.12S.ncbi.${DATE}.ids_novel.txt > target_species.12S.ncbi.${DATE}.novel.fa
        echo "Rename"
        seqkit replace -p '(.+)$' -r '{kv}' -k target_species.12S.ncbi.${DATE}.rename.txt target_species.12S.ncbi.${DATE}.novel.fa --keep-key > target_species.12S.ncbi.${DATE}.novel_relabel.fa
        echo "Concatenate"
        cat ${label}.raw.fa target_species.12S.ncbi.${DATE}.novel_relabel.fa > ${label}.plus_target.fa
        mv ${label}.plus_target.fa ${label}.raw.fa
      else
        echo "  No gap-filling requested, moving on..."
      fi
      # pcr - allow ~10% mismatch per primer
      usearch -search_pcr ${label}.raw.fa -db ${SCRIPTS}/12S_primers.fa -strand both -maxdiffs 4 -minamp 420 -maxamp 470 -ampout ${label}.amp.fa
      if [ -s ${label}.amp.fa ]
      then
        # remove primers
        usearch -fastx_truncate ${label}.amp.fa -stripleft 30 -stripright 28 -fastaout ${label}.amp_only.fa
        # blast
        makeblastdb -in ${label}.amp_only.fa -dbtype nucl
        blastn -db ${label}.amp_only.fa -query ${label}.raw.fa -out ${label}.amp.blastn -task blastn -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -num_threads 2
        cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=360){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords
        # remove blastdb
        rm ${label}.amp_only.fa.n*
        # get amplicons
        seqtk subseq ${label}.raw.fa ${label}.amp.blastn.coords | awk -F: '{print $1}' > ${label}.amp_blast.fa
        # check for ambigs
        seqtk comp ${label}.amp_blast.fa | awk '{sum=0; for(i=7;i<=9; i++) {sum+=$i} {perc=100*sum/$2} if(perc==0.00) print $1}' - | seqtk subseq ${label}.amp_blast.fa - > ${label}.amp_blast.noN.fa
        # align
        mafft --adjustdirection --retree 2 --reorder --thread 4 ${label}.amp_blast.noN.fa > ${label}.amp_blast.noN.mafft.fa
      else
        echo "  No sequences with matching primer sequences found, skipping..."
        echo ""
        continue
      fi
      #cleanup
      mkdir -p ./intermediate_files
      mv ${label}.keep_relabel.fa ./intermediate_files
      mv ${label}.raw.fa ./intermediate_files
      mv ${label}.amp_blast.fa ./intermediate_files
      mv ${label}.amp_blast.noN.fa ./intermediate_files
      mv ${label}.amp_blast.noN.mafft.fa ./intermediate_files
      rm ${label}*
      rm -f target_species.*
    elif [ ${locus} = "lrRNA" ]
    then
      echo "  Amplicon should be ${locus}..."
      # add extras if requested
      if [ ${EXTRAS} = "yes" ] && [ -s extra_16S.fa ]
      then
        echo "  Adding extra sequences..."
        # add extra sequences (make sure on single line and uppercase)
        cat ${label}.keep_relabel.fa <(awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' extra_16S.fa) | seqtk seq -l0 - > ${label}.raw.fa
      else
        echo "  No extras to add. Skipping..."
        cp ${label}.keep_relabel.fa ${label}.raw.fa
      fi
      if [ -s ${GAPFILL} ]
      then
        echo "  Checking for target species..."
        splist=($(cut -f1 ${GAPFILL}))
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"l-rRNA\" OR \"16S rRNA\" OR \"16S ribosomal RNA\"\)" | efilter -location mitochondrion -source genbank | efetch -format fasta >> target_species.16S.genbank.${DATE}.fa; done`
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"l-rRNA\" OR \"16S rRNA\" OR \"16S ribosomal RNA\"\)" | efilter -location mitochondrion -source refseq | efetch -format fasta >> target_species.16S.refseq.${DATE}.fa; done`
        cat target_species.16S.genbank.${DATE}.fa target_species.16S.refseq.${DATE}.fa > target_species.16S.ncbi.${DATE}.fa
        # edit NC_ to NC (and AC_ to AC)
        sed -i 's/\(>[A-Z][A-Z]\)\(_\)/\1/g' target_species.16S.ncbi.${DATE}.fa
        #discard dupicates based on ID
        awk '/^>/{f=!d[$1];d[$1]=1}f' target_species.16S.ncbi.${DATE}.fa > tmp
        #drop UNVERIFIED and any with {}s in header
        seqtk subseq tmp <(grep ">" tmp | sed 's/>//g' | sed '/UNVERIFIED/d' | sed '/{/d' | sed '/}/d') > target_species.16S.ncbi.${DATE}.fa
        rm tmp
        # get current headers and make a renaming file to convert to same format as protax database (Genus_species_accession)
        grep ">" target_species.16S.ncbi.${DATE}.fa | sed 's/>//g' > target_species.16S.ncbi.${DATE}.headers.txt
        cut -f1,2,3 -d" " target_species.16S.ncbi.${DATE}.headers.txt | sed 's/\.[0-9]/ /g' | awk '{print $1 "\t" $2 "_" $3 "_" $1}' > target_species.16S.ncbi.${DATE}.rename.txt
        # simplify current headers
        awk '{print $1}' target_species.16S.ncbi.${DATE}.fa | sed 's/\(^>[A-Z0-9].\+\)\(\.[0-9]\)/\1/g' > target_species.16S.ncbi.${DATE}.ids.fa
        # get list of accessions and check against list of previously seen accessions
        grep ">" target_species.16S.ncbi.${DATE}.ids.fa | sed 's/>//g' | sort -u | cut -f1 -d" " | cut -f1 -d"." > target_species.16S.ncbi.${DATE}.ids
        ids=($(cut -f1 target_species.16S.ncbi.${DATE}.ids))
        for id in ${ids[@]}
        do
          if grep -q ${id} ${label}.raw.fa
          then
            printf '%s\t%s\n' ${id} "true" >> target_species.16S.ncbi.${DATE}.ids_status.txt
          else
            printf '%s\t%s\n' ${id} "false" >> target_species.16S.ncbi.${DATE}.ids_status.txt
          fi
        done
        # select previously unseen sequences
        grep "false" target_species.16S.ncbi.${DATE}.ids_status.txt | cut -f1 > target_species.16S.ncbi.${DATE}.ids_novel.txt
        seqtk subseq target_species.16S.ncbi.${DATE}.ids.fa target_species.16S.ncbi.${DATE}.ids_novel.txt > target_species.16S.ncbi.${DATE}.novel.fa
        echo "Rename"
        seqkit replace -p '(.+)$' -r '{kv}' -k target_species.16S.ncbi.${DATE}.rename.txt target_species.16S.ncbi.${DATE}.novel.fa --keep-key > target_species.16S.ncbi.${DATE}.novel_relabel.fa
        echo "Concatenate"
        cat ${label}.raw.fa target_species.16S.ncbi.${DATE}.novel_relabel.fa > ${label}.plus_target.fa
        mv ${label}.plus_target.fa ${label}.raw.fa
      else
        echo "  No gap-filling requested, moving on..."
      fi
      # pcr - allow ~10% mismatch per primer
      usearch -search_pcr ${label}.raw.fa -db ${SCRIPTS}/16S_primers.fa -strand both -maxdiffs 2 -minamp 130 -maxamp 140 -ampout ${label}.amp.fa
      if [ -s ${label}.amp.fa ]
      then
        # remove primers
        usearch -fastx_truncate ${label}.amp.fa -stripleft 19 -stripright 21 -fastaout ${label}.amp_only.fa
        # blast
        makeblastdb -in ${label}.amp_only.fa -dbtype nucl
        blastn -db ${label}.amp_only.fa -query ${label}.raw.fa -out ${label}.amp.blastn -task blastn -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -num_threads 2
        cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=90){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords
        # remove blastdb
        rm ${label}.amp_only.fa.n*
        # get amplicons
        seqtk subseq ${label}.raw.fa ${label}.amp.blastn.coords | awk -F: '{print $1}' > ${label}.amp_blast.fa
        # check for ambigs
        seqtk comp ${label}.amp_blast.fa | awk '{sum=0; for(i=7;i<=9; i++) {sum+=$i} {perc=100*sum/$2} if(perc==0.00) print $1}' - | seqtk subseq ${label}.amp_blast.fa - > ${label}.amp_blast.noN.fa
        # align
        mafft --adjustdirection --retree 2 --reorder --thread 4 ${label}.amp_blast.noN.fa > ${label}.amp_blast.noN.mafft.fa
      else
        echo "  No sequences with matching primer sequences found, skipping..."
        echo ""
        continue
      fi
      #cleanup
      mkdir -p ./intermediate_files
      mv ${label}.keep_relabel.fa ./intermediate_files
      mv ${label}.raw.fa ./intermediate_files
      mv ${label}.amp_blast.fa ./intermediate_files
      mv ${label}.amp_blast.noN.fa ./intermediate_files
      mv ${label}.amp_blast.noN.mafft.fa ./intermediate_files
      rm ${label}*
      rm -f target_species.*
    elif [ ${locus} = "Cytb" ]
    then
      echo "  Amplicon should be ${locus}..."
      # add extras if requested
      if [ ${EXTRAS} = "yes" ] && [ -s extra_Cytb.fa ]
      then
        echo "  Adding extra sequences..."
        # add extra sequences (make sure on single line and uppercase)
        cat ${label}.keep_relabel.fa <(awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' extra_Cytb.fa) | seqtk seq -l0 - > ${label}.raw.fa
      else
        echo "  No extras to add. Skipping..."
        cp ${label}.keep_relabel.fa ${label}.raw.fa
      fi
      if [ -s ${GAPFILL} ]
      then
        echo "  Checking for target species..."
        splist=($(cut -f1 ${GAPFILL}))
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"cytb\" OR \"CytB\" OR \"cytochrome b\"\)" | efilter -location mitochondrion -source genbank | efetch -format fasta >> target_species.Cytb.genbank.${DATE}.fa; done`
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"cytb\" OR \"CytB\" OR \"cytochrome b\"\)" | efilter -location mitochondrion -source refseq | efetch -format fasta >> target_species.Cytb.refseq.${DATE}.fa; done`
        cat target_species.Cytb.genbank.${DATE}.fa target_species.Cytb.refseq.${DATE}.fa > target_species.Cytb.ncbi.${DATE}.fa
        # edit NC_ to NC (and AC_ to AC)
        sed -i 's/\(>[A-Z][A-Z]\)\(_\)/\1/g' target_species.Cytb.ncbi.${DATE}.fa
        #discard dupicates based on ID
        awk '/^>/{f=!d[$1];d[$1]=1}f' target_species.Cytb.ncbi.${DATE}.fa > tmp
        #drop UNVERIFIED and any with {}s in header
        seqtk subseq tmp <(grep ">" tmp | sed 's/>//g' | sed '/UNVERIFIED/d' | sed '/{/d' | sed '/}/d') > target_species.Cytb.ncbi.${DATE}.fa
        rm tmp
        # get current headers and make a renaming file to convert to same format as protax database (Genus_species_accession)
        grep ">" target_species.Cytb.ncbi.${DATE}.fa | sed 's/>//g' > target_species.Cytb.ncbi.${DATE}.headers.txt
        cut -f1,2,3 -d" " target_species.Cytb.ncbi.${DATE}.headers.txt | sed 's/\.[0-9]/ /g' | awk '{print $1 "\t" $2 "_" $3 "_" $1}' > target_species.Cytb.ncbi.${DATE}.rename.txt
        # simplify current headers
        awk '{print $1}' target_species.Cytb.ncbi.${DATE}.fa | sed 's/\(^>[A-Z0-9].\+\)\(\.[0-9]\)/\1/g' > target_species.Cytb.ncbi.${DATE}.ids.fa
        # get list of accessions and check against list of previously seen accessions
        grep ">" target_species.Cytb.ncbi.${DATE}.ids.fa | sed 's/>//g' | sort -u | cut -f1 -d" " | cut -f1 -d"." > target_species.Cytb.ncbi.${DATE}.ids
        ids=($(cut -f1 target_species.Cytb.ncbi.${DATE}.ids))
        for id in ${ids[@]}
        do
          if grep -q ${id} ${label}.raw.fa
          then
            printf '%s\t%s\n' ${id} "true" >> target_species.Cytb.ncbi.${DATE}.ids_status.txt
          else
            printf '%s\t%s\n' ${id} "false" >> target_species.Cytb.ncbi.${DATE}.ids_status.txt
          fi
        done
        # select previously unseen sequences
        grep "false" target_species.Cytb.ncbi.${DATE}.ids_status.txt | cut -f1 > target_species.Cytb.ncbi.${DATE}.ids_novel.txt
        seqtk subseq target_species.Cytb.ncbi.${DATE}.ids.fa target_species.Cytb.ncbi.${DATE}.ids_novel.txt > target_species.Cytb.ncbi.${DATE}.novel.fa
        echo "Rename"
        seqkit replace -p '(.+)$' -r '{kv}' -k target_species.Cytb.ncbi.${DATE}.rename.txt target_species.Cytb.ncbi.${DATE}.novel.fa --keep-key > target_species.Cytb.ncbi.${DATE}.novel_relabel.fa
        echo "Concatenate"
        cat ${label}.raw.fa target_species.Cytb.ncbi.${DATE}.novel_relabel.fa > ${label}.plus_target.fa
        mv ${label}.plus_target.fa ${label}.raw.fa
      else
        echo "  No gap-filling requested, moving on..."
      fi
      # pcr - allow ~10% mismatch per primer
      usearch -search_pcr ${label}.raw.fa -db ${SCRIPTS}/Cytb_primers.fa -strand both -maxdiffs 5 -minamp 365 -maxamp 385 -ampout ${label}.amp.fa
      if [ -s ${label}.amp.fa ]
      then
        # remove primers
        usearch -fastx_truncate ${label}.amp.fa -stripleft 30 -stripright 31 -fastaout ${label}.amp_only.fa
        # blast
        makeblastdb -in ${label}.amp_only.fa -dbtype nucl
        blastn -db ${label}.amp_only.fa -query ${label}.raw.fa -out ${label}.amp.blastn -task blastn -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -num_threads 2
        cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=300){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords
        # remove blastdb
        rm ${label}.amp_only.fa.n*
        # get amplicons
        seqtk subseq ${label}.raw.fa ${label}.amp.blastn.coords | awk -F: '{print $1}' > ${label}.amp_blast.fa
        # check for ambigs
        seqtk comp ${label}.amp_blast.fa | awk '{sum=0; for(i=7;i<=9; i++) {sum+=$i} {perc=100*sum/$2} if(perc==0.00) print $1}' - | seqtk subseq ${label}.amp_blast.fa - > ${label}.amp_blast.noN.fa
        # align
        mafft --adjustdirection --retree 2 --reorder --thread 4 ${label}.amp_blast.noN.fa > ${label}.amp_blast.noN.mafft.fa
      else
        echo "  No sequences with matching primer sequences found, skipping..."
        echo ""
        continue
      fi
      #cleanup
      mkdir -p ./intermediate_files
      mv ${label}.keep_relabel.fa ./intermediate_files
      mv ${label}.raw.fa ./intermediate_files
      mv ${label}.amp_blast.fa ./intermediate_files
      mv ${label}.amp_blast.noN.fa ./intermediate_files
      mv ${label}.amp_blast.noN.mafft.fa ./intermediate_files
      rm ${label}*
      rm -f target_species.*
    elif [ ${locus} = "COI" ]
    then
      echo "  Amplicon should be ${locus}..."
      # add extras if requested
      if [ ${EXTRAS} = "yes" ] && [ -s extra_COI.fa ]
      then
        echo "  Adding extra sequences..."
        # add extra sequences (make sure on single line and uppercase)
        cat ${label}.keep_relabel.fa <(awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' extra_COI.fa) | seqtk seq -l0 - > ${label}.raw.fa
      else
        echo "  No extras to add. Skipping..."
        cp ${label}.keep_relabel.fa ${label}.raw.fa
      fi
      if [ -s ${GAPFILL} ]
      then
        echo "  Checking for target species..."
        splist=($(cut -f1 ${GAPFILL}))
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"COX1\" OR \"COI\" OR \"cytochrome oxidase c subunit I\"\)" | efilter -location mitochondrion -source genbank | efetch -format fasta >> target_species.COI.genbank.${DATE}.fa; done`
        `for sp in ${splist[@]}; do esearch -db nucleotide -query "${sp} [ORGN] AND \(\"COX1\" OR \"COI\" OR \"cytochrome oxidase c subunit I\"\)" | efilter -location mitochondrion -source refseq | efetch -format fasta >> target_species.COI.refseq.${DATE}.fa; done`
        cat target_species.COI.genbank.${DATE}.fa target_species.COI.refseq.${DATE}.fa > target_species.COI.ncbi.${DATE}.fa
        # edit NC_ to NC (and AC_ to AC)
        sed -i 's/\(>[A-Z][A-Z]\)\(_\)/\1/g' target_species.COI.ncbi.${DATE}.fa
        #discard dupicates based on ID
        awk '/^>/{f=!d[$1];d[$1]=1}f' target_species.COI.ncbi.${DATE}.fa > tmp
        #drop UNVERIFIED and any with {}s in header
        seqtk subseq tmp <(grep ">" tmp | sed 's/>//g' | sed '/UNVERIFIED/d' | sed '/{/d' | sed '/}/d') > target_species.COI.ncbi.${DATE}.fa
        rm tmp
        # get current headers and make a renaming file to convert to same format as protax database (Genus_species_accession)
        grep ">" target_species.COI.ncbi.${DATE}.fa | sed 's/>//g' > target_species.COI.ncbi.${DATE}.headers.txt
        cut -f1,2,3 -d" " target_species.COI.ncbi.${DATE}.headers.txt | sed 's/\.[0-9]/ /g' | awk '{print $1 "\t" $2 "_" $3 "_" $1}' > target_species.COI.ncbi.${DATE}.rename.txt
        # simplify current headers
        awk '{print $1}' target_species.COI.ncbi.${DATE}.fa | sed 's/\(^>[A-Z0-9].\+\)\(\.[0-9]\)/\1/g' > target_species.COI.ncbi.${DATE}.ids.fa
        # get list of accessions and check against list of previously seen accessions
        grep ">" target_species.COI.ncbi.${DATE}.ids.fa | sed 's/>//g' | sort -u | cut -f1 -d" " | cut -f1 -d"." > target_species.COI.ncbi.${DATE}.ids
        ids=($(cut -f1 target_species.COI.ncbi.${DATE}.ids))
        for id in ${ids[@]}
        do
          if grep -q ${id} ${label}.raw.fa
          then
            printf '%s\t%s\n' ${id} "true" >> target_species.COI.ncbi.${DATE}.ids_status.txt
          else
            printf '%s\t%s\n' ${id} "false" >> target_species.COI.ncbi.${DATE}.ids_status.txt
          fi
        done
        # select previously unseen sequences
        grep "false" target_species.COI.ncbi.${DATE}.ids_status.txt | cut -f1 > target_species.COI.ncbi.${DATE}.ids_novel.txt
        seqtk subseq target_species.COI.ncbi.${DATE}.ids.fa target_species.COI.ncbi.${DATE}.ids_novel.txt > target_species.COI.ncbi.${DATE}.novel.fa
        echo "Rename"
        seqkit replace -p '(.+)$' -r '{kv}' -k target_species.COI.ncbi.${DATE}.rename.txt target_species.COI.ncbi.${DATE}.novel.fa --keep-key > target_species.COI.ncbi.${DATE}.novel_relabel.fa
        echo "Concatenate"
        cat ${label}.raw.fa target_species.COI.ncbi.${DATE}.novel_relabel.fa > ${label}.plus_target.fa
        mv ${label}.plus_target.fa ${label}.raw.fa
      else
        echo "  No gap-filling requested, moving on..."
      fi
      # pcr - allow ~10% mismatch per primer
      usearch -search_pcr ${label}.raw.fa -db ${SCRIPTS}/COI_primers.fa -strand both -maxdiffs 3 -minamp 240 -maxamp 260 -ampout ${label}.amp.fa
      if [ -s ${label}.amp.fa ]
      then
        # remove primers
        usearch -fastx_truncate ${label}.amp.fa -stripleft 26 -stripright 19 -fastaout ${label}.amp_only.fa
        # blast
        makeblastdb -in ${label}.amp_only.fa -dbtype nucl
        blastn -db ${label}.amp_only.fa -query ${label}.raw.fa -out ${label}.amp.blastn -task blastn -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -num_threads 2
        cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=195){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords
        # remove blastdb
        rm ${label}.amp_only.fa.n*
        # get amplicons
        seqtk subseq ${label}.raw.fa ${label}.amp.blastn.coords | awk -F: '{print $1}' > ${label}.amp_blast.fa
        # check for ambigs
        seqtk comp ${label}.amp_blast.fa | awk '{sum=0; for(i=7;i<=9; i++) {sum+=$i} {perc=100*sum/$2} if(perc==0.00) print $1}' - | seqtk subseq ${label}.amp_blast.fa - > ${label}.amp_blast.noN.fa
        # align
        mafft --adjustdirection --retree 2 --reorder --thread 4 ${label}.amp_blast.noN.fa > ${label}.amp_blast.noN.mafft.fa
      else
        echo "  No sequences with matching primer sequences found, skipping..."
        echo ""
        continue
      fi
      #cleanup
      mkdir -p ./intermediate_files
      mv ${label}.keep_relabel.fa ./intermediate_files
      mv ${label}.raw.fa ./intermediate_files
      mv ${label}.amp_blast.fa ./intermediate_files
      mv ${label}.amp_blast.noN.fa ./intermediate_files
      mv ${label}.amp_blast.noN.mafft.fa ./intermediate_files
      rm ${label}*
      rm -f target_species.*
    else
      echo "Locus not recognised."
      exit 1
    fi
  done
  #cleanup any leftovers from skipped loci
  rm -f MIDORI_*.fa
  rm -f MIDORI_*.txt
  rm -f target_species*
  #cleanup input files
  mv MIDORI_UNIQUE_* ./intermediate_files
  [ -f ${GAPFILL} ] && mv ${GAPFILL} ./intermediate_files
  [ -f extra_12S.fa ] && mv extra_12S.fa ./intermediate_files
  [ -f extra_16S.fa ] && mv extra_16S.fa ./intermediate_files
  [ -f extra_Cytb.fa ] && mv extra_Cytb.fa ./intermediate_files
  [ -f extra_COI.fa ] && mv extra_COI.fa ./intermediate_files
  echo ""
  end=`date +%s`
  runtime=$((end-start))
  echo ""
  echo "Module 1 took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
  echo ""
  echo "Module 1 complete. Stopping now for manual inspection of alignments *.mafft.fa inside ./intermediate_files."
  echo "Restart script when happy with alignments (save as *.mafft_edit.fa in present directory even if no edits are made)."
  echo "Input files have been moved to ./intermediate_files"
  echo ""
  echo "Enjoy :-)"
  echo ""
  }

  function module_two {
    echo "Starting with manually checked alignments..."
    echo ""
    echo "  Looking for alignments..."
    if ls MIDORI_*.mafft_edit.fa 1> /dev/null 2>&1
    then
    FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
    echo "  Done."
    else
      echo "  No matching FASTAs found."
      exit 1
    fi
    #get working taxonomy file (original from get_taxonomy.sh)
    tr '\r' '\n' < ${TAXON}_ITIS_taxonomy.txt > tmp
    mv tmp ${TAXON}_ITIS_taxonomy.txt
    sed 's/ /_/g' -i ${TAXON}_ITIS_taxonomy.txt
    #edit _R_ labels from mafft if any
    for file in ${FILES[@]}
    do
      sed -i 's/>_R_/>/g' ${file}
    done
    #get full list of species & accessions
    cat MIDORI_*.mafft_edit.fa | grep '>' | sed 's/>//g' | sort -u | sed 's/_/\t/g' | awk '{print $1 "_" $2 "_" $3 "\t" $1 "_" $2}' > MIDORI.raw_id2sp.txt
    cat MIDORI_*.mafft_edit.fa | grep '>' | sed 's/>//g' | sort -u | sed 's/_/\t/g' | awk '{print $1 "_" $2 "_" $3 "\t" $3}' > MIDORI.raw_id2acc.txt
    # get list pf uniq species
    awk '{print $2}' MIDORI.raw_id2sp.txt | sort -u > MIDORI.raw_sp.txt
    # non-matching
    tabtk isct -1 7 -2 1 -c ${TAXON}_ITIS_taxonomy.txt  MIDORI.raw_sp.txt | sed 's/_/ /g' > MIDORI_${TAXON}.ITIS_mismatch_sp.txt
    # remove misbehavers
    sed -i '/Hemidactylus adensis/d' MIDORI_${TAXON}.ITIS_mismatch_sp.txt
    sed -i '/Hemidactylus awashensis/d' MIDORI_${TAXON}.ITIS_mismatch_sp.txt
    # lookup non-matching in CoL in first instance, lookup failures in ITIS, output list of failures to discard
    Rscript --vanilla ${SCRIPTS}/get_taxonomy_mismatches.R ${TAXON}
    # combine CoL and ITIS outputs (if ITIS exists)
    if [ -s ${TAXON}_ITIS_matched_taxonomy.txt ]
    then
      cat ${TAXON}_CoL_matched_taxonomy.txt ${TAXON}_ITIS_matched_taxonomy.txt > ${TAXON}_matched_taxonomy.txt
      tr '\r' '\n' < ${TAXON}_matched_taxonomy.txt > tmp
      sed 's/ /_/g' tmp > ${TAXON}_matched_taxonomy.txt
      rm tmp
    else
      tr '\r' '\n' < ${TAXON}_CoL_matched_taxonomy.txt > tmp
      sed 's/ /_/g' tmp > ${TAXON}_matched_taxonomy.txt
      rm tmp
    fi
    # double check there are only target taxa present
    if [ ${TAXON} = "Tetrapoda" ]
    then
      grep -E "Amphibia|Aves|Mammalia|Reptilia" ${TAXON}_matched_taxonomy.txt > tmp
      mv tmp ${TAXON}_matched_taxonomy.txt
    else
      grep ${TAXON} ${TAXON}_matched_taxonomy.txt > tmp
      mv tmp ${TAXON}_matched_taxonomy.txt
    fi
    cat <(awk '{print $0 "\t" "ITIS_sp" "\t" $7}' ${TAXON}_ITIS_taxonomy.txt) <(awk '{print $0}' ${TAXON}_matched_taxonomy.txt) | sed '1ikingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstatus\tquery' | sed 's/ /_/g' > ${TAXON}.combined_taxonomy.txt
    # update missing list
    sed -i 's/ /_/g' ${TAXON}.ITIS_missing_sp.txt
    tabtk isct -1 7 -2 1 -c ${TAXON}.combined_taxonomy.txt <(tabtk isct -1 9 -2 1 -c ${TAXON}.combined_taxonomy.txt <(cat MIDORI.raw_sp.txt ${TAXON}.ITIS_missing_sp.txt)) > MIDORI_${TAXON}.missing_sp.txt
    if [ -s MIDORI_${TAXON}.missing_sp.txt ]
    then
      sed -i '$a\Hemidactylus_adensis\' MIDORI_${TAXON}.missing_sp.txt
      sed -i '$a\Hemidactylus_awashensis\' MIDORI_${TAXON}.missing_sp.txt
      # check if missing sp have genus in taxonomy, if so add a new line, if not check if known synonym and add a new line, if not flag for deletion
      splist=($(cut -f1 MIDORI_${TAXON}.missing_sp.txt))
      for sp in ${splist[@]}
      do
        genus=($(echo ${sp} | cut -f1 -d'_'))
        if grep -q ${genus} <(cut -f6 ${TAXON}.combined_taxonomy.txt)
        then
          printf '%s\t%s\t%s\n' ${sp} "present" "genus_match" >> MIDORI.genus_status.txt
          tabtk grep -f7 "${genus}_" ${TAXON}.combined_taxonomy.txt | cut -f1-6 | sort -u | awk -v sp="${sp}" '{print $0 "\t" sp "\t" "genus_match" "\t" sp}' >> genus_match_taxo
        else
          if grep -q ${genus} <(cut -f9 ${TAXON}.combined_taxonomy.txt)
          then
            printf '%s\t%s\t%s\n' ${sp} "absent" "synonym_match" >> MIDORI.genus_status.txt
            species=($(echo ${sp} | cut -f2 -d'_'))
            tabtk grep -f9 "${genus}_" ${TAXON}.combined_taxonomy.txt | cut -f1-6 | sort -u | awk -v species="${species}" -v sp="${sp}" '{print $0 "\t" $6 "_" species "\t" "genus_synonym_match" "\t" sp}' >> genus_match_taxo
          else
            printf '%s\t%s\t%s\n' ${sp} "absent" "no_match" >> MIDORI.genus_status.txt
            printf '%s\n' ${sp} >> ${TAXON}.missing_sp_to_delete.txt
          fi
        fi
      done
      # add newly created genus-match species taxonomy to existing
      cat ${TAXON}.combined_taxonomy.txt genus_match_taxo > tmp
      mv tmp ${TAXON}.combined_taxonomy.txt
      #cleanup
      mv ${TAXON}_ITIS_taxonomy.txt ./intermediate_files
      rm -f ${TAXON}.ITIS_missing_sp.txt
      rm -f MIDORI.genus_status.txt
      rm -f MIDORI_${TAXON}.missing_sp.txt
      rm -f genus_match_taxo
      rm -f ${TAXON}_CoL_matched_taxonomy.txt
      rm -f ${TAXON}_ITIS_matched_taxonomy.txt
      rm -f ${TAXON}_matched_taxonomy.txt
      rm -f MIDORI.raw_sp.txt
      rm -f MIDORI_${TAXON}.ITIS_mismatch_sp.txt
    else
      #cleanup
      mv ${TAXON}_ITIS_taxonomy.txt ./intermediate_files
      rm -f ${TAXON}.ITIS_missing_sp.txt
      rm -f MIDORI_${TAXON}.missing_sp.txt
      rm -f ${TAXON}_CoL_matched_taxonomy.txt
      rm -f ${TAXON}_ITIS_matched_taxonomy.txt
      rm -f ${TAXON}_matched_taxonomy.txt
      rm -f MIDORI.raw_sp.txt
      rm -f MIDORI_${TAXON}.ITIS_mismatch_sp.txt
      echo "No failed lookups to check, going straight to Sativa step."
      echo ""
      module_three
    fi
    echo ""
    end=`date +%s`
    runtime=$((end-start))
    echo ""
    echo "Module 2 took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
    echo ""
    echo "Module 2 complete. Stopping now for manual inspection of failed species lookups (in ${TAXON}.missing_sp_to_delete.txt)."
    echo "If a failed lookup can be resolved, remove from ${TAXON}.missing_sp_to_delete.txt and add taxonomy to a tab-delimited file named ${TAXON}.missing_sp_taxonomy.txt with columns for kingdom,phylum,class,order,family,genus,species,status,query - 'status' should be something short and descriptive ("_" instead of spaces; eg. "mispelling" or "manual_synonym") and 'query' should be the entry in ${TAXON}.missing_sp_to_delete.txt. ${TAXON}.missing_sp_taxonomy.txt must not have a header line when the script is restarted."
    echo "If all failed lookups are resolved, delete ${TAXON}.missing_sp_to_delete.txt. If some/all failed lookups cannot be resolved, keep the relevant species names in ${TAXON}.missing_sp_to_delete.txt. When restarting the script it will check for the presence of this file and act accordingly (sequences for these species will be discarded)."
    echo "If no failed lookups can be resolved, do no create ${TAXON}.missing_sp_taxonomy.txt, leave ${TAXON}.missing_sp_to_delete.txt as it is."
    echo "Restart script when happy."
    echo ""
    echo "Enjoy :-)"
    echo ""
  }
  function module_three {
    echo "Merging manually checked taxonomy with existing where relevant and discarding failed lookups..."
    echo ""
    if [ -s ${TAXON}.missing_sp_taxonomy.txt ] && ! [ -s ${TAXON}.missing_sp_to_delete.txt ]
    then
      echo "  ${TAXON}.missing_sp_taxonomy.txt present, ${TAXON}.missing_sp_to_delete.txt absent..."
      cat ${TAXON}.combined_taxonomy.txt ${TAXON}.missing_sp_taxonomy.txt > ${TAXON}.consensus_taxonomy.txt
      # edit alignment names
      FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
      for file in ${FILES[@]}
      do
        label=($(echo ${file} | sed 's/.amp_blast.noN.mafft_edit.fa//g'))
        cp ${file} ${label}.final.fa
      done
    elif [ -s ${TAXON}.missing_sp_taxonomy.txt ] && [ -s ${TAXON}.missing_sp_to_delete.txt ]
    then
      echo "  Both ${TAXON}.missing_sp_taxonomy.txt and ${TAXON}.missing_sp_to_delete.txt present..."
      cat ${TAXON}.combined_taxonomy.txt ${TAXON}.missing_sp_taxonomy.txt > ${TAXON}.consensus_taxonomy.txt
      # discard failed lookups
      sed -i 's/ /_/g' ${TAXON}.missing_sp_to_delete.txt
      ids=($(cut -f1 MIDORI.raw_id2acc.txt))
      for i in ${ids[@]}
      do
        sp=($(echo ${i} | cut -f1,2 -d'_'))
        if grep -q ${sp} ${TAXON}.missing_sp_to_delete.txt
        then
          printf '%s\t%s\n' ${i} "delete" >> MIDORI.taxonomy_status.txt
        else
          printf '%s\t%s\n' ${i} "keep" >> MIDORI.taxonomy_status.txt
          printf '%s\n' ${i} >> MIDORI.seqs_to_retain.txt
        fi
      done
      #need to delete seqs by selecting ones to keep
      FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
      for file in ${FILES[@]}
      do
        label=($(echo ${file} | sed 's/.amp_blast.noN.mafft_edit.fa//g'))
        seqtk subseq ${file} MIDORI.seqs_to_retain.txt > ${label}.final.fa
      done
    elif ! [ -s ${TAXON}.missing_sp_taxonomy.txt ] && [ -s ${TAXON}.missing_sp_to_delete.txt ]
    then
      echo "  ${TAXON}.missing_sp_taxonomy.txt absent, ${TAXON}.missing_sp_to_delete.txt present..."
      cp ${TAXON}.combined_taxonomy.txt ${TAXON}.consensus_taxonomy.txt
      # discard failed lookups
      sed -i 's/ /_/g' ${TAXON}.missing_sp_to_delete.txt
      ids=($(cut -f1 MIDORI.raw_id2acc.txt))
      for i in ${ids[@]}
      do
        sp=($(echo ${i} | cut -f1,2 -d'_'))
        if grep -q ${sp} ${TAXON}.missing_sp_to_delete.txt
        then
          printf '%s\t%s\n' ${i} "delete" >> MIDORI.taxonomy_status.txt
        else
          printf '%s\t%s\n' ${i} "keep" >> MIDORI.taxonomy_status.txt
          printf '%s\n' ${i} >> MIDORI.seqs_to_retain.txt
        fi
      done
      #need to delete seqs by selecting ones to keep
      FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
      for file in ${FILES[@]}
      do
        label=($(echo ${file} | sed 's/.amp_blast.noN.mafft_edit.fa//g'))
        seqtk subseq ${file} MIDORI.seqs_to_retain.txt > ${label}.final.fa
      done
    else
      echo "  No ${TAXON}.missing_sp_taxonomy.txt present, assume this is becuase there were no failed lookups to check, moving on..."
      cp ${TAXON}.combined_taxonomy.txt ${TAXON}.consensus_taxonomy.txt
      # edit alignment names
      FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
      for file in ${FILES[@]}
      do
        label=($(echo ${file} | sed 's/.amp_blast.noN.mafft_edit.fa//g'))
        cp ${file} ${label}.final.fa
      done
    fi
    echo "Getting final species names labels..."
    # get protax & sativa names
    join -1 2 -2 9 -o 1.1,1.3,2.7 <(cut -f1 MIDORI.raw_id2sp.txt | awk 'BEGIN{FS="_"}{print $1 "_" $2 "_" $3 "\t" $1 "_" $2 "\t" $3}' | sort -k2,2) <(sort -k9,9 ${TAXON}.consensus_taxonomy.txt) > MIDORI.final_id2acc2sp.txt
    awk '{print $1 "\t" $2}' MIDORI.final_id2acc2sp.txt > ${TAXON}.final_rename_seqs_sativa.txt
    echo "Generating master sativa taxonomy file..."
    join -1 3 -2 7 -o 1.2,2.2,2.3,2.4,2.5,2.6,2.7 <(sort -k3,3 MIDORI.final_id2acc2sp.txt) <(sort -k7,7 ${TAXON}.consensus_taxonomy.txt) | awk '{print $1 "\t" "Eukaryota;" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7}' | sed 's/_/ /g' | sort -u > ${TAXON}.final_taxonomy_sativa.txt
    echo "Working on each locus seperately..."
    FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.final.fa" | sed 's/\.\///g'))
    for file in ${FILES[@]}
    do
      label=($(echo ${file} | sed 's/.final.fa//g'))
      locus=($(echo ${label} | sed 's/*_//g'))
      echo "  Working on ${locus}..."
      echo "  Renaming alignment and generating sativa taxonomy file..."
      #rename alignments for sativa
      seqkit replace -p '(.+)$' -r '{kv}' -k ${TAXON}.final_rename_seqs_sativa.txt ${file} --keep-key > ${label}.final_for_sativa.fa
      #make sativa taxonomy
      join -1 1 -2 1 <(grep ">" ${label}.final_for_sativa.fa | sed 's/>//g' | sort -k1,1) <(sort -k1,1 ${TAXON}.final_taxonomy_sativa.txt) | awk '{print $1 "\t" $2 " " $3}' > ${label}.final_for_sativa.tax
    done
    # sativa
    echo "  Running sativa for each alignment..."
    FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.final_for_sativa.fa" | sed 's/\.\///g'))
    for file in ${FILES[@]}
    do
      label=($(echo ${file} | sed 's/.final_for_sativa.fa//g'))
      locus=($(echo ${label} | sed 's/*_//g'))
      mkdir ./${label}_sativa
      sativa -s ${file} -t ${label}.final_for_sativa.tax -x ZOO -T 4 -n ${label} -o ./${label}_sativa
      #cleanup
      mv ${label}.final.fa ./intermediate_files
      rm ${label}.final_for_sativa.tax
    done
    #cleanup
    mv MIDORI_*.amp_blast.noN.mafft_edit.fa ./intermediate_files
    rm ${TAXON}.final_rename_seqs_sativa.txt
    rm MIDORI.taxonomy_status.txt
    rm MIDORI.seqs_to_retain.txt
    rm MIDORI.final_id2acc2sp.txt
    rm -f ${TAXON}.missing_sp_taxonomy.txt
    rm ${TAXON}.missing_sp_to_delete.txt
    rm ${TAXON}.consensus_taxonomy.txt
    rm ${TAXON}.combined_taxonomy.txt
    rm MIDORI.raw_id2sp.txt
    rm MIDORI.raw_id2acc.txt
    echo ""
    end=`date +%s`
    runtime=$((end-start))
    echo ""
    echo "Module 3 took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
    echo ""
    echo "Module 3 complete. Stopping now for manual inspection of mislabelled sequences in ./MIDORI_locus_sativa/MIDORI_locus.mis"
    echo "To skip manual editing, do nothing and restart script."
    echo "To make changes, create file ./MIDORI_locus_sativa/MIDORI_locus.mis_to_delete and copy all confirmed mislabelled sequence accessions as a single column, tab-delimited list. These will be deleted at the start of module 4."
    echo "For sequences where species-level or genus-level mislabelling can be resolved, make corrections directly in ${TAXON}.final_taxonomy_sativa.txt (i.e. replace the taxonomic classification for that sequence with the correct one), this will be used to rename sequences."
    echo "Make higher level changes to the taxonomy at your own risk - untested."
    echo "Restart script when happy."
    echo ""
    echo "Enjoy :-)"
    echo ""
  }
  function module_four {
    FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.final_for_sativa.fa" | sed 's/\.\///g'))
    for file in ${FILES[@]}
    do
      label=($(echo ${file} | sed 's/.final_for_sativa.fa//g'))
      if [ -s ./${label}_sativa/${label}.mis_to_delete ]
      then
        # select accepted
        cut -f1 ./${label}_sativa/${label}.mis_to_delete | sort -k1,1 > ${label}.sativa_flagged.txt
        tabtk isct -c -1 1 -1 1 ${label}.sativa_flagged.txt <(grep ">" ${file} | sed 's/>//g') | seqtk subseq ${file} - > ${label}.final_clean.fa
        # rename
        sed 's/;/\t/g' ${TAXON}.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > ${TAXON}.final_rename_seqs_protax.txt
        seqkit replace -p '(.+)$' -r '{kv}' -k ${TAXON}.final_rename_seqs_protax.txt ${label}.final_clean.fa --keep-key > ${label}.final_clean_relabel.fa
        seqbuddy ${label}.final_clean_relabel.fa --clean_seq > ${label}.final_clean_relabel.unalign.fa
        # make final taxonomy
        cut -f2 ${TAXON}.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > ${TAXON}.final_protax_taxonomy.txt
        # get list of sp
        grep ">" ${label}.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > ${label}.sp_list.txt
        # for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
        touch ${label}.uniq.fa
        splist=($(cut -f1 ${label}.sp_list.txt))
        for sp in ${splist[@]}
        do
          grep ${sp} ${label}.final_clean_relabel.unalign.fa | sed 's/>//g' | seqtk subseq ${label}.final_clean_relabel.unalign.fa - > ${sp}.fa
          count=($(grep -c ">" ${sp}.fa))
          if [ $count -gt 1 ]
          then
            mafft --thread 2 --retree 2 --reorder ${sp}.fa > ${sp}.mafft.fa
            perl ${SCRIPTS}/collapsetypes_v4.6.pl -infile=${sp}.mafft.fa -nrdiffs=0
            seqbuddy ${sp}.mafft.fa.unique_haplotypes --clean_seq > ${sp}.uniq.fa
            cat ${label}.uniq.fa ${sp}.uniq.fa > tmp
          else
            cat ${label}.uniq.fa ${sp}.fa > tmp
          fi
          mv tmp ${label}.uniq.fa
          rm ${sp}*
        done
        #cleanup
        mv ${label}.final_clean_relabel.unalign.fa ./intermediate_files
        rm ${label}.sativa_flagged.txt
        rm ${label}.final_clean.fa
        rm ${label}.final_clean_relabel.fa
        rm ${label}.sp_list.txt
      else
        # select accepted
        sed -e '1d;2d;3d;4d;5d' ./${label}_sativa/${label}.mis | cut -f1 | sort -k1,1 > ${label}.sativa_flagged.txt
        tabtk isct -c -1 1 -1 1 ${label}.sativa_flagged.txt <(grep ">" ${file} | sed 's/>//g') | seqtk subseq ${file} - > ${label}.final_clean.fa
        # rename
        sed 's/;/\t/g' ${TAXON}.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > ${TAXON}.final_rename_seqs_protax.txt
        seqkit replace -p '(.+)$' -r '{kv}' -k ${TAXON}.final_rename_seqs_protax.txt ${label}.final_clean.fa --keep-key > ${label}.final_clean_relabel.fa
        seqbuddy ${label}.final_clean_relabel.fa --clean_seq > ${label}.final_clean_relabel.unalign.fa
        # make final taxonomy
        cut -f2 ${TAXON}.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > ${TAXON}.final_protax_taxonomy.txt
        # get list of sp
        grep ">" ${label}.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > ${label}.sp_list.txt
        # for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
        touch ${label}.uniq.fa
        splist=($(cut -f1 ${label}.sp_list.txt))
        for sp in ${splist[@]}
        do
          grep ${sp} ${label}.final_clean_relabel.unalign.fa | sed 's/>//g' | seqtk subseq ${label}.final_clean_relabel.unalign.fa - > ${sp}.fa
          count=($(grep -c ">" ${sp}.fa))
          if [ $count -gt 1 ]
          then
            mafft --thread 2 --retree 2 --reorder ${sp}.fa > ${sp}.mafft.fa
            perl ${SCRIPTS}/collapsetypes_v4.6.pl -infile=${sp}.mafft.fa -nrdiffs=0
            seqbuddy ${sp}.mafft.fa.unique_haplotypes --clean_seq > ${sp}.uniq.fa
            cat ${label}.uniq.fa ${sp}.uniq.fa > tmp
          else
            cat ${label}.uniq.fa ${sp}.fa > tmp
          fi
          mv tmp ${label}.uniq.fa
          rm ${sp}*
        done
        #cleanup
        mv ${label}.final_for_sativa.fa ./intermediate_files
        mv ${label}.final_clean_relabel.unalign.fa ./intermediate_files
        rm ${label}.sativa_flagged.txt
        rm ${label}.final_clean.fa
        rm ${label}.final_clean_relabel.fa
        rm ${label}.sp_list.txt
      fi
      # end up with mix of upper/lower case and one/multi-line seqs, change to single-line uppercase
      awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' ${label}.uniq.fa | seqtk seq -l0 - > tmp
      locus=$(echo $label | cut -f2 -d "_")
      if [ ${locus} = "srRNA" ]
      then
        locus="12S"
      elif [ ${locus} = "lrRNA" ]
      then
        locus="16S"
      fi
      mv tmp ${TAXON}.final_database.${locus}.fa
    done
    echo ""
    #cleanup
    mv ${TAXON}.final_taxonomy_sativa.txt ./intermediate_files
    mv ./MIDORI_*_sativa/ ./intermediate_files
    rm ${TAXON}.final_rename_seqs_protax.txt
    rm MIDORI_*.uniq.fa
    #end
    end=`date +%s`
    runtime=$((end-start))
    echo ""
    echo "Module 4 took `echo $runtime | awk '{printf "%.2f", $1/3600}'` hours"
    echo ""
    echo "Module 4 complete. You have reached the end of get_sequences.sh"
    echo ""
    echo "Final database sequences are in ${TAXON}.final_database.locus.fa, final taxonomy file is in ${TAXON}.final_protax_taxonomy.txt"
    echo ""
    echo "Next step: train PROTAX models with either:"
    echo "  - train_protax.sh for unweighted models"
    echo "  - train_weighted_protax.sh for models weighted using a list of expected species"
  	echo ""
  	echo "Have a nice day :-)"
    echo ""
  }

##### MAIN
  if [ ${MODULE} = "one" ]
  then
    echo "Starting at Module 1"
    start=`date +%s`
    echo "The time now is $(date)."
    echo "Information from you:"
    echo "Local FASTA files to be added: ${EXTRAS}"
    echo "Targeted gap-filling: ${GAPFILL}"
    echo "Target taxon: ${TAXON}"

    echo ""
    module_one
  elif [ ${MODULE} = "two" ]
  then
    echo "Restarting at Module 2"
    start=`date +%s`
    echo "The time now is $(date)."
    echo ""
    module_two
  elif [ ${MODULE} = "three" ]
  then
    echo "Restarting at Module 3"
    start=`date +%s`
    echo "The time now is $(date)."
    echo ""
    module_three
  elif [ ${MODULE} = "four" ]
  then
    echo "Restarting at Module 4"
    start=`date +%s`
    echo "The time now is $(date)."
    echo ""
    module_four
  else
    echo "Value of module not recognised - must be 'one', 'two', 'three' or 'four'."
    exit 1
  fi
