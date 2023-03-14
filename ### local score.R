### local score
## Detecting genomic outliers in a sequence of correlated scores, from Fariello et al 2017.
## Detecting minor and major QTL in GWAS using the local score approach, from Bonhomme et al. 2019.
setwd("/net/cephfs/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

library(tidyverse)

### Import libraries.
library(ggplot2)
library(data.table)
library(RColorBrewer)

### import data
phen=read.table("pheno_T_H_B.txt", h=T,na.strings=".")
phen <- phen[-6]

mrk=read.table("./subsetT_H_B_geno.bim") # a table with first column the chromosomes, and the second column the Position of each marker0 The order should be the same than .map file, seems like we just want col 1 and 4
mrk <- mrk %>%
        select(c(1,4)) 

colnames(mrk)=c("chrom", "po")

traits=colnames(phen)[2:ncol(phen)]

for (trait in traits){

#trait <- "Heightd34"
#windowsFonts(
#A=windowsFont("Calibri")
#)
chrom   POS     pval    log_pval        cor     lindley th      signif
Local_Score_NbFlowers_xi2.txt

data=read.table(paste("./T_H_B_",trait,".emmax.ps",sep=""), h=T)
data <- data %>%
        select(-1) 
colnames(data)=c("SNP", "Beta", "pvalue")

emmax <- bind_cols(mrk$chrom, mrk$po, data$Beta, data$pvalue)  

emmax <- emmax %>% 
        dplyr::rename(
        chrom = ...1,
        POS =  ...2,
        Beta =  ...3,
        pval =  ...4) 



emmax <- emmax %>% 
        mutate(ID = paste(emmax$chrom, emmax$POS, sep = '_'))

emmax$manhattan=0

emmax$pval <- as.numeric(unlist(emmax$pval))
emmax <- as.data.table(emmax)
setkey(emmax, chrom)
emmax$pval[emmax$pval==0]=1e-16#min(mydata$pval[-which(mydata$pval ==0)])
Nchr=length(emmax[,unique(chrom)])
emmax$log_pval = -log10(emmax$pval)
#Import functions.
source('/net/cephfs/shares/schiestl.systbot.uzh/variants/data/tfig/scorelocalfunctions.R')

### Computation of absolute position in the genome. 
#This is useful for doing genomewide plots. 
chrInfo=emmax[,.(L=.N,cor=autocor(pval)),chrom]
setkey(chrInfo,chrom)
data.table(chr=emmax[,unique(chrom),], S=cumsum(c(0,chrInfo$L[-Nchr])))

#    chr       S
# 1:   1       0
# 2:   2  502582
# 3:   3 1035712
# 4:   4 1645638
# 5:   5 2017172
# 6:   6 2502878
# 7:   7 3034952
# 8:   8 3549049
# 9:   9 3957440
#10:  10 4741068
### Choice of $\xi$ (1,2,3 or 4)
#emmax <- cbind(emmax, score=NA, lindley=NA)


# The score mean must be negative, ksi must be chosen between mean and max of -log10(p-value)
# xi should be between 2 and 3 for gwas, make sure that mean score is negative
xi=2
emmax[,score:= emmax$log_pval-xi]
mean(emmax$score)
emmax[,lindley:=lindley(score),chrom]

mean(emmax$log_pval);max(emmax$log_pval)

#hist(mydata$score)
max(emmax$lindley)

# Compute significance threshold for each chromosome
## Uniform distribution of p-values
chrInfo[,th:=thresUnif(L, cor, xi),chrom]
(emmax=emmax[chrInfo])


(sigZones=emmax[chrInfo,sig_sl(lindley, POS, unique(th)),chrom])
sigZones$QTL_length <- sigZones$end - sigZones$beg
sigZones <- filter(sigZones,peak>0)
sigZonesSNPDens = data.frame()
all_SNP_signif_log4 = vector()

emmax$signif = round(emmax$lindley-emmax$th,3)
emmax$th = round(emmax$th,3)
emmax$cor = round(emmax$cor,3)
emmax$log_pval = round(emmax$log_pval,6)


  for (zone in c(1:nrow(sigZones))) {
   CHR = sigZones$chr[zone]
   BEG = sigZones$beg[zone]
   END = sigZones$end[zone]
   SNP_number_log4 = nrow(emmax[emmax$chr==CHR & emmax$ps>=BEG & emmax$ps<=END & emmax$log_plrt>4])
   SNP_signif_log4 = emmax$rs[emmax$chr==CHR & emmax$ps>=BEG & emmax$ps<=END & emmax$log_plrt>4]
   SNP_number_log5 = nrow(emmax[emmax$chr==CHR & emmax$ps>=BEG & emmax$ps<=END & emmax$log_plrt>5])
   SNP_number_log6 = nrow(emmax[emmax$chr==CHR & emmax$ps>=BEG & emmax$ps<=END & emmax$log_plrt>6])
   subsetzone = cbind(sigZones[zone,],SNP_number_log4,SNP_number_log5,SNP_number_log6)
   sigZonesSNPDens = rbind(sigZonesSNPDens,subsetzone)
   all_SNP_signif_log4 = append(all_SNP_signif_log4,SNP_signif_log4)
  }

### write results into table
results = as.data.frame(cbind(emmax$chrom,emmax$POS,emmax$pval,emmax$log_pval,emmax$cor,emmax$lindley,emmax$th,emmax$signif))
colnames(results) = c("chrom","POS","pval","log_pval","cor","lindley","th","signif")

### so here just recall emmax = results
#emmax <- results

## can just plot right from results to not load in again
write.table(results,file=paste0("./Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F,)
write.table(as.data.frame(sigZonesSNPDens),file=paste0("SigZones_Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F)
}

emmax$manhattan=0

### create plotting data frames
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


### plot it: ylab = -log10(pval) with sigZones in red
tiff(paste("./T_NH_B_Manhattan_local_score_", trait,"_xi", xi, "_sigZones.tiff",sep=""), res=600, width=25, height=15, unit="cm",compression="lzw") 

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
  plot.window(xlim=c(deb_emmax,fin_emmax),ylim=c(0,max(max(emmax$lindley,max(emmax$th))*1.05)))
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
tiff(paste("./T_NH_B_Manhattan_local_score_", trait,"_xi", xi, "_lindley.tiff",sep=""), res=600, width=25, height=15, unit="cm",compression="lzw") 

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
