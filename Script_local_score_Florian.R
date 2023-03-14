## Detecting genomic outliers in a sequence of correlated scores, from Fariello et al 2017.
## Detecting minor and major QTL in GWAS using the local score approach, from Bonhomme et al. 2019.

### Import libraries.
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(data.table)
library(tidyverse)


#Import functions
# download functions in https://forge-dga.jouy.inra.fr/documents/809

### run the command without any changes from here to "####end function"
source('./scorelocalfunctions.R')

sig_zone_local_score <- function(lindley,pos,th){
  zones=c(0,0,0,0)	
  LIND=lindley
  auxpos=pos
  while(max(LIND,na.rm=TRUE)>=th){
    M_loc <- which.max(LIND)
    if(min(LIND[1:M_loc],na.rm=TRUE)==0){
      beg_peak <- max(which(LIND[1:M_loc]==min(LIND[1:M_loc],na.rm=TRUE)))+1 #found first 0 before max peak to identify peak start
    }else{
      beg_peak <- max(which(LIND[1:M_loc]==min(LIND[1:M_loc],na.rm=TRUE)))
    }
    if(length(which(LIND[M_loc+1:length(LIND)]==0))>1){
      end_peak <- min(which(LIND[M_loc+1:length(LIND)]==min(LIND[M_loc+1:length(LIND)],na.rm=TRUE)))+M_loc #found first 0 after max peak to identify peak end
    }else{
      end_peak <- M_loc
    }
    zones=rbind(zones, c(auxpos[beg_peak],auxpos[end_peak],auxpos[M_loc],max(LIND)))
    LIND <- LIND[-(beg_peak:end_peak)]
    auxpos <- auxpos[-(beg_peak:end_peak)]
  }
  zones=matrix(zones, ncol=4)
  zones=data.table(beg=zones[,1],end=zones[,2],max_zone=zones[,3],peak=zones[,4])
  if (nrow(zones)>1){zones=zones[-1,]}
  return(zones)
}


files = list.files("./output/", pattern=".assoc.txt", all.files=FALSE, full.names=FALSE)
files = tools::file_path_sans_ext(files)
files = tools::file_path_sans_ext(files)

####end function

# start script for your analyses
mrk=read.table("./list_chrom_pos_order_fam.txt", h=T) # the one you generated during GWA (I don't remember how you call it)

# for the loop, write in vector traits your different variables
#traits=c("....")
#for (trait in traits[1:length(traits)]){
  
trait="NbFlowers" 
  ### Choice of $\xi$ (1,2,3 or 4)
  #xi=2
  xi=3
  
  data=read.table(paste("./resistance_",trait,".emmax.ps",sep=""), h=T)
  colnames(data)=c("SNP", "Beta", "pvalue")
  data=cbind(mrk$chrom, mrk$pos, data$Beta, data$pvalue)
  colnames(data)=c("chr", "pos", "Beta", "pvalue")
  mydata=as.data.table(data[,c(which(colnames(data)=="chr"), which(colnames(data)=="pos"), which(colnames(data)=="pval"))])
 
  setkey(mydata, chr) # for sorting mydata according to chr, the column will be marked as sorted with a key
  mydata$pval[mydata$pval==0]=1e-16
  Nchr=length(mydata[,unique(chr)])
    
  ### Computation of absolute position in the genome. 
    #This is useful for doing genomewide plots. 
    chrInfo=mydata[,.(L=.N,cor=autocor(pval)),chr]
    setkey(chrInfo,chr)
    data.table(chr=mydata[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr])))
    mydata$log_pval = -log10(mydata$pval)
    
    mydata[,score:= mydata$log_pval-xi]
    mean(mydata$score)
    mydata[,lindley:=lindley(score),chr]
    
    # The score mean must be negative, ksi must be chosen between mean and max of -log10(p-value)
    mean(mydata$log_pval);max(mydata$log_pval)
    mean(mydata$score)
    #hist(mydata$score)
    max(mydata$lindley)
    
    # Compute significance threshold for each chromosome
    ## Uniform distribution of p-values
    chrInfo[,th:=thresUnif(L, cor, xi),chr]
    mydata=mydata[chrInfo]
    sigZones <- mydata[chrInfo,sig_zone_local_score(lindley, pos, unique(th)),chr]
    sigZones$QTL_length <- sigZones$end - sigZones$beg
    
    
    sigZones <- filter(sigZones,peak>0)
    
    sigZonesSNPDens = data.frame()
    all_SNP_signif_log4 = vector()
    
    mydata$signif = round(mydata$lindley-mydata$th,3)
    mydata$th = round(mydata$th,3)
    mydata$cor = round(mydata$cor,3)
    mydata$log_pval = round(mydata$log_pval,6)
    
    for (zone in c(1:nrow(sigZones))) {
     CHR = sigZones$chr[zone]
     BEG = sigZones$beg[zone]
     END = sigZones$end[zone]
     SNP_number_log4 = nrow(mydata[mydata$chr==CHR & mydata$ps>=BEG & mydata$ps<=END & mydata$log_plrt>4])
     SNP_signif_log4 = mydata$rs[mydata$chr==CHR & mydata$ps>=BEG & mydata$ps<=END & mydata$log_plrt>4]
     SNP_number_log5 = nrow(mydata[mydata$chr==CHR & mydata$ps>=BEG & mydata$ps<=END & mydata$log_plrt>5])
     SNP_number_log6 = nrow(mydata[mydata$chr==CHR & mydata$ps>=BEG & mydata$ps<=END & mydata$log_plrt>6])
     subsetzone = cbind(sigZones[zone,],SNP_number_log4,SNP_number_log5,SNP_number_log6)
     sigZonesSNPDens = rbind(sigZonesSNPDens,subsetzone)
     all_SNP_signif_log4 = append(all_SNP_signif_log4,SNP_signif_log4)
    }
    
    
## write your results    
    results = as.data.frame(cbind(mydata$chr,mydata$pos,mydata$BFdB,mydata$log_pval,mydata$cor,mydata$lindley,mydata$th,mydata$signif))
    colnames(results) = c("chr","pos","pval", "log_pval","cor","lindley","th","signif")
    
    write.table(results,file=paste0("./Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F,)
    write.table(as.data.frame(sigZonesSNPDens),file=paste0("SigZones_Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F)

    #    }  # to run for loop only


  
############
### Plot it 
####
trait="NbFlowers" 
xi=2  

emmax=read.table(paste("./Local_Score_",trait,"_xi", xi, ".txt",sep=""), h=T)
emmax$manhattan=0

for (c in 1:length(unique(emmax$chrom))){
  assign(paste("emmax",c,sep=""),emmax[emmax[,1]==unique(emmax$chrom)[c],c(1:ncol(emmax))])
}

emmax1$manhattan=emmax1$POS
emmax2$manhattan=emmax2$POS+max(emmax1$manhattan)
emmax3$manhattan=emmax3$POS+max(emmax2$manhattan)
emmax4$manhattan=emmax4$POS+max(emmax3$manhattan)
emmax5$manhattan=emmax5$POS+max(emmax4$manhattan)
emmax6$manhattan=emmax6$POS+max(emmax5$manhattan)
emmax7$manhattan=emmax7$POS+max(emmax6$manhattan)
emmax8$manhattan=emmax8$POS+max(emmax7$manhattan)
emmax9$manhattan=emmax9$POS+max(emmax8$manhattan)
emmax10$manhattan=emmax10$POS+max(emmax9$manhattan)

emmax1$col="grey20"
emmax2$col="grey40"
emmax3$col="grey20"
emmax4$col="grey40"
emmax5$col="grey20"
emmax6$col="grey40"
emmax7$col="grey20"
emmax8$col="grey40"
emmax9$col="grey20"
emmax10$col="grey40"

emmax$manhattan=c(emmax1$manhattan,emmax2$manhattan,emmax3$manhattan,emmax4$manhattan,emmax5$manhattan,emmax6$manhattan,emmax7$manhattan,emmax8$manhattan,emmax9$manhattan,emmax10$manhattan)
emmax$col=c(emmax1$col,emmax2$col,emmax3$col,emmax4$col,emmax5$col,emmax6$col,emmax7$col,emmax8$col,emmax9$col,emmax10$col)


deb_emmax=min(emmax$manhattan)
fin_emmax=max(emmax$manhattan)


#windowsFonts(
#  A=windowsFont("Calibri")
)


### plot it: ylab = -log10(pval) with sigZones in red
tiff(paste("./Manhattan_local_score_", trait,"_xi", xi, "_sigZones.tiff",sep=""), res=600, width=25, height=15, unit="cm",compression="lzw") 

m=matrix(1:3,3,1)
layout(m,c(5),c(1.5,10,3))
layout.show(3)

par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
  text(0.05,0.2,paste(trait, sep=""),cex=1.5)
par(mar=c(0,2,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  plot(emmax$manhattan,-log10(emmax$pval),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb_emmax,fin_emmax),axes=FALSE, cex =0.3 ,  col=emmax$col ,pch=16)
  SNP_signif_log4 = filter(emmax, emmax$signif>0, emmax$log_pval>4)
  points(SNP_signif_log4$manhattan,-log10(SNP_signif_log4$pval), col="red", cex=1, pch=16)
  axis(2,cex.axis=1,lwd=0.5,tck=0.03)
  plot.window(xlim=c(deb_emmax,fin_emmax),ylim=c(0,max(max(emmax$lind,max(emmax$th))*1.05)))
  lines(emmax$manhattan,emmax$lindley,pch=16,col="red",lwd=1)
  for(chrom in unique(emmax$chrom)){
  segments(x0=min(emmax$manhattan[emmax$chrom==chrom]) ,y0= unique(emmax$th[emmax$chrom==chrom]),
           x1=max(emmax$manhattan[emmax$chrom==chrom]),y1=unique(emmax$th[emmax$chrom==chrom]),
           lty=3, col="red")
}  
box(lwd=0.5)
axis(4,col="red",lwd=0.5,tck=0.03)
mtext("-log10(pval)", side=2, line=1, cex=0.8)

dev.off()



### plot it with upper panel ylab = -log10(pval), below panel ylab = lindley score


tiff(paste("./Manhattan_local_score_", trait,"_xi", xi, "_lindley.tiff",sep=""), res=600, width=25, height=15, unit="cm",compression="lzw") 

m=matrix(1:5,5,1)
layout(m,c(5),c(1.5,10,1.5,10,3))
layout.show(5)

par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
  text(0.05,0.2,paste(trait, sep=""),cex=1.5)
par(mar=c(0,2,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  plot(emmax$manhattan, -log10(emmax$pval),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb_emmax,fin_emmax),axes=FALSE, cex =0.3 ,  col=emmax$col ,pch=16)
  axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb_emmax,fin_emmax,10),labels=F, tck=0)
  axis(2,cex.axis=1,lwd=0.5,tck=0.03)
  box(lwd=0.5)
  mtext("pval", side=2, line=1, cex=0.8)
par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
par(mar=c(0,2,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  plot(emmax$manhattan,emmax$lindley,xlab="",ylab="",frame.plot=FALSE,xlim=c(deb_emmax,fin_emmax),axes=FALSE, cex =0.5 ,  col=emmax$col ,pch=16)
  axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb_emmax,fin_emmax,10),labels=F, tck=0)
  axis(2,cex.axis=1,lwd=0.5,tck=0.03)
  #abline(h=mean(emmax$th), col="grey", lty=2, las=1)
  box(lwd=0.5)
  mtext("lindley score", side=2, line=1, cex=0.8)
par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)

dev.off()

