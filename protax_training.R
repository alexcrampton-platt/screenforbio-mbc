#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("Welcome to protax_training.R")
message("")
message("Step 1: Load necessary packages etc")
message("")
protaxdir<-args[1]
taxon<-args[2]
loci<-list.files(pattern=paste0(taxon,".final_database"))
for(i in 1:length(loci)){
  loci[i]<-sub(".fa","",loci[i])
  loci[i]<-sub(paste0(taxon,".final_database."),"",loci[i])
}
source(paste0(protaxdir,"/amcmc.rcode.txt"))
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)
num.params=1+4
ind=1001:2000
message("")
message("Step 2: Run four iterations of training for each marker")
message("This will take some time...")
message("")
for(locus in loci){
  folder=paste0("./model_",locus,"/")
  message(paste0("Working on ",locus," in folder ",folder))
  for(level in c(1,2,3,4)){
    message(paste0("Working on level",level))
    dat=read.xdata(paste0(folder,"train",level,".xdat"))
    ppa=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
    initstate=initialize.adaptation(ppa$params[2000,])
    ppb=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
    initstate=initialize.adaptation(ppb$params[2000,])
    ppc=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
    initstate=initialize.adaptation(ppc$params[2000,])
    ppd=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
    pdf(paste0(folder,"training_plot_",locus,"_level",level,"a_MCMC.pdf"))
    traceplot.all(ppa,ind,num.levels=1, title="iter1")
    amcmc.diagnostic.plot(ppa)
    dev.off()
    pdf(paste0(folder,"training_plot_",locus,"_level",level,"b_MCMC.pdf"))
    traceplot.all(ppb,ind,num.levels=1, title="iter2")
    amcmc.diagnostic.plot(ppb)
    dev.off()
    pdf(paste0(folder,"training_plot_",locus,"_level",level,"c_MCMC.pdf"))
    traceplot.all(ppc,ind,num.levels=1, title="iter3")
    amcmc.diagnostic.plot(ppc)
    dev.off()
    pdf(paste0(folder,"training_plot_",locus,"_level",level,"d_MCMC.pdf"))
    traceplot.all(ppd,ind,num.levels=1, title="iter4")
    amcmc.diagnostic.plot(ppd)
    dev.off()
    k=which.max(ppa$postli[ind])
    write.postparams(ppa,paste0(folder,"mcmc",level,"a"),ind[k])
    k=which.max(ppb$postli[ind])
    write.postparams(ppb,paste0(folder,"mcmc",level,"b"),ind[k])
    k=which.max(ppc$postli[ind])
    write.postparams(ppc,paste0(folder,"mcmc",level,"c"),ind[k])
    k=which.max(ppd$postli[ind])
    write.postparams(ppd,paste0(folder,"mcmc",level,"d"),ind[k])
  }
}
message("")
message("End of protax_training.R")
q(save="no")
