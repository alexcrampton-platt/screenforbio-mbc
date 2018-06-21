#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("Welcome to get_taxonomy.R")
message("")
message("Step 1: Check for and load library(taxize)")
message("")
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
message("Step 2: Get list of species with downstream()")
taxon<-args[1]
rank<-args[2]
taxon_list<-downstream(get_tsn(paste(taxon),rank=paste(rank),kingdom="Animalia"),db="itis",downto="species")
taxon_sp<-as.data.frame(taxon_list[1])[,1]
message("Step 3: Get taxonomy with classification()")
taxon_class<-classification(as.tsn(taxon_sp),db="itis",return_id=FALSE)
message("")
message("Step 4: Write lookups to file")
taxon_tab<-cbind(taxon_class)
write.table(cbind(rep("Eukaryota",nrow(taxon_tab)),taxon_tab$phylum,taxon_tab$class,taxon_tab$order,taxon_tab$family,taxon_tab$genus,taxon_tab$species),paste0(taxon,"_ITIS_taxonomy.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
message("	File is:")
message("")
paste0(taxon,"_ITIS_taxonomy.txt")
message("")
message("End of get_taxonomy.R")
q(save="no")
