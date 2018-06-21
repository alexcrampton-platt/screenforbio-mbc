args = commandArgs(trailingOnly=TRUE)

message("Welcome to training_plots.R")
message("")
model<-args[1]
taxon<-args[2]
locus<-args[3]
modtype<-args[4]
setwd(paste0(model))
message("Step 1: Make plots")
aplot = function(prob,correct,add=F,col="black",name="") {
  n=length(prob)
  if (!add) {
   plot(0,type="n",xlab="cumulative prob %",ylab="cumulative correct %",xlim=c(0,100),ylim=c(0,100),las=1)
   grid()
   abline(0,1,col="gray")
   title(sprintf("%s\n%d items, correct %.1f %%",name,n,sum(correct)/n*100))
  }
  s=sort(prob,dec=F,index.return=T)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  lines(x,y,col=col)
  points(x[n],y[n],pch=19,cex=1,col=col)
}
message(" Done.")
message("Step 2: Save plots")
nimi=c("order","family","genus","species")
pdf(paste0(modtype,"_",taxon,"_",locus,"_biasaccuracy.pdf"),width=8,height=8)
par(mfrow=c(2,2))
for (i in 1:4) {
 file=sprintf("query%d.cor",i)
 a=read.table(file,header=F)
 aplot(a[,2],a[,3],name=sprintf("%s: %s",taxon,nimi[i]))
}
dev.off()
message(" Done.")
message("")
message("End of training_plots.R")
quit(save="no")
