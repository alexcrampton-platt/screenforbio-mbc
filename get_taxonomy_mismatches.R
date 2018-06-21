#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("Welcome to get_taxonomy_mismatches.R")
message("")
message("Step1: Check for and load library(taxize)")
pkgLoad <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE, repos='https://cloud.r-project.org/')
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
pkgLoad("taxize")
message("")
message("Step 2: Import list of species")
taxon<-args[1]
mismatch<-read.table(paste0("MIDORI_",taxon,".ITIS_mismatch_sp.txt"),header=F,sep="\t")$V1
message("")
message("Step 3: Get CoL matches with classification()")
mismatch_class<-classification(mismatch,rank="species",db="col",return_id=FALSE,rows=1)
message("")
message("Step 4: Write successful CoL lookups to file")
col_matched<-mismatch_class[unlist(lapply(mismatch_class,function(x)is.character(x[[1]])))]
class(col_matched) <- "classification"
col_tab<-cbind(col_matched)
write.table(cbind(rep("Eukaryota",nrow(col_tab)),col_tab$phylum,col_tab$class,col_tab$order,col_tab$family,col_tab$genus,col_tab$species,rep("CoL_match",nrow(col_tab)),col_tab$query),paste0(taxon,"_CoL_matched_taxonomy.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
message("	Done. File is:")
paste0(taxon,"_CoL_matched_taxonomy.txt")
message("")
message("Step 5: Get ITIS matches with classification()")
col_unmatched <- mismatch_class[unlist(lapply(mismatch_class,function(x)is.logical(x[[1]])))]
col_missing <- names(col_unmatched)
missing_class<-classification(col_missing,rank="species",db="itis",return_id=FALSE,rows=1)
message("")
message("Step 6: Write failed lookups to file")
itis_unmatched <- missing_class[unlist(lapply(missing_class,function(x)is.logical(x[[1]])))]
itis_missing <- names(itis_unmatched)
if((length(itis_missing))!=0){
write.table(itis_missing,paste0(taxon,".ITIS_missing_sp.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	message("	File is:")
	paste0(taxon,".ITIS_missing_sp.txt")
} else {
	message("	No failed lookups to write out.")
}
message("")
message("Step 7: Write successful ITIS lookups to file")
itis_matched<-missing_class[unlist(lapply(missing_class,function(x)is.character(x[[1]])))]
class(itis_matched) <- "classification"
if (length(itis_matched)!=0){
    itis_tab<-cbind(itis_matched)
    write.table(cbind(rep("Eukaryota",nrow(itis_tab)),itis_tab$phylum,itis_tab$class,itis_tab$order,itis_tab$family,itis_tab$genus,itis_tab$species,rep("ITIS_match",nrow(itis_tab)),itis_tab$query),paste0(taxon,"_ITIS_matched_taxonomy.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    message("	Done. File is:")
    paste0(taxon,"_ITIS_matched_taxonomy.txt")
    message("")
}else{
   message("No matches found in ITIS")
}
message("End of get_taxonomy_mismatches.R")
q(save="no")
