### combine chrom files into one and take local score for drift sims, plot or pca after

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

setwd("/shares/schiestl.systbot.uzh/variants/data/tfig")

if (!require(learnPopGen)) install.packages('learnPopGen')
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

chrome <- c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10")
n<- chrome[7]
subset = "total"

# faster to just load in with data table, rather than convert to data.table object
# single read in, just to try
res <- fread(paste0(subset, "_chrome", n, "_1_fdr.tsv")) 

# if need to combine from multiple files
res10 <- fread("total_chromeA10_1_fdr.tsv", select = c(1,2,10005))
res_list <- list(res, res10)
res <- rbindlist(res_list) 

## use unix to merge and then load
## awk goes line by line so maybe faster elsewhere with cat or grep
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print $1, $2, $10005}
' subset_AB_T_B_chromeA0*_1_fdr.tsv subset_AB_T_B_chromeA10_1_fdr.tsv> subset_AB_T_B_chrome_onlysig_fdr.tsv

# if done with awk, then has headers inside, remove and check
res <- fread("total_drift_pval_fin.tsv")

res$CHROM <- as.factor(res$CHROM)
levels(res$CHROM)

setkey(res, CHROM)
res <- res[!"CHROM"]
res$CHROM<- droplevels(res$CHROM)
levels(res$CHROM)

res$POS <- as.numeric(res$POS)

# see if things are ordered before
keycol <-c("CHROM","POS")
res<- setorderv(res, keycol)


### do local score with drift sim, then can plot sig zones for pca
res$pval <- as.numeric(unlist(res$pval))

setkey(res, CHROM)
Nchr=length(res[,unique(CHROM)])

### so here, the pvalue is probably not 0 -16, it may create false peaks, 
#res$pvalue[res$pvalue==0]=1e-16

min(res$pval[-which(res$pval ==0)])
#[1] 1e-04 - so lets set it to - 0
res$pvalue[res$pvalue==0]=1e-04

res$log_pval = -log10(res$pvalue)

#Import functions.
source('/shares/schiestl.systbot.uzh/variants/data/tfig/scorelocalfunctions.R')

chrInfo=res[,.(L=.N,cor=autocor(pval)),CHROM]
setkey(chrInfo,CHROM)
data.table(chr=res[,unique(CHROM),], S=cumsum(c(0,chrInfo$L[-Nchr])))

#    chr       S
# 1:   4       0


# The score mean must be negative, ksi must be chosen between mean and max of -log10(p-value)
# xi should be between 2 and 3 for gwas, make sure that mean score is negative
xi=2
res[,score:= res$log_pval-xi]
mean(res$score)
res[,lindley:=lindley(score),CHROM]

mean(res$log_pval);max(res$log_pval)

#hist(mydata$score)
max(res$lindley)

# Compute significance threshold for each chromosome
## Uniform distribution of p-values
chrInfo[,th:=thresUnif(L, cor, xi),CHROM]
(res=res[chrInfo])


(sigZones=res[chrInfo,sig_sl(lindley, as.numeric(POS), unique(th)),CHROM])
sigZones$QTL_length <- sigZones$end - sigZones$beg
sigZones <- filter(sigZones,peak>0)
sigZonesSNPDens = data.frame()
all_SNP_signif_log4 = vector()

res$signif = round(res$lindley-res$th,3)
res$th = round(res$th,3)
res$cor = round(res$cor,3)
res$log_pval = round(res$log_pval,6)


for (zone in c(1:nrow(sigZones))) {
 CHR = sigZones$chr[zone]
 BEG = sigZones$beg[zone]
 END = sigZones$end[zone]
 SNP_number_log4 = nrow(res[res$chr==CHR & res$ps>=BEG & res$ps<=END & res$log_plrt>4])
 SNP_signif_log4 = res$rs[res$chr==CHR & res$ps>=BEG & res$ps<=END & res$log_plrt>4]
 SNP_number_log5 = nrow(res[res$chr==CHR & res$ps>=BEG & res$ps<=END & res$log_plrt>5])
 SNP_number_log6 = nrow(res[res$chr==CHR & res$ps>=BEG & res$ps<=END & res$log_plrt>6])
 subsetzone = cbind(sigZones[zone,],SNP_number_log4,SNP_number_log5,SNP_number_log6)
 sigZonesSNPDens = rbind(sigZonesSNPDens,subsetzone)
 all_SNP_signif_log4 = append(all_SNP_signif_log4,SNP_signif_log4)
}

### write results into table, leave as emmax and other cols so easier to plot
emmax = as.data.frame(cbind(res$CHROM,res$POS,res$pvalue,res$log_pval,res$cor,res$lindley,res$th,res$signif))
colnames(emmax) = c("chrom","POS","pval","log_pval","cor","lindley","th","signif")

# not numeric for some reason so convert 
emmax$POS <- as.numeric(unlist(emmax$POS))
emmax$pval <- as.numeric(unlist(emmax$pval))
emmax$log_pval <- as.numeric(unlist(emmax$log_pval))
emmax$cor <- as.numeric(unlist(emmax$cor))
emmax$lindley <- as.numeric(unlist(emmax$lindley))
emmax$th <- as.numeric(unlist(emmax$th))
emmax$signif <- as.numeric(unlist(emmax$signif))

## can just plot right from results to not load in again
#write.table(emmax,file=paste0("./Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F,)
#write.table(as.data.frame(sigZonesSNPDens),file=paste0("SigZones_Local_Score_",trait,"_xi", xi, ".txt"),quote=F, sep="\t",row.names=F)


### so here just recall emmax = results
#emmax <- results

### plot the local score
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
  text(0.05,0.2,paste(subset, sep=""),cex=1.5)
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
  text(0.05,0.2,paste(subset, sep=""),cex=1.5)
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


### below is only if we use e1-16, we went with 1 -04 which makes more sense, do not need to remove >0 because non
### then we export sigzones and do vcftools on it
ext2<- sigZones[,c(1:3)]
fwrite(ext2, sep = "\t", paste0(subset, "_subset_1e04_snps_vcf.bed"))
## after use script_vcftools_extract
## then run the script_plink_pca_sig_sub



# lets say that we extract all sig zones where qtl is bigger than 0
ext <- sigZones[QTL_length > 0]

## check to see that the beg and end = pos
#sub[POS == 4269204 & CHROM == "A07"] # it does so lets extract within the zone

ext$CHROM <- droplevels(ext$CHROM)
levels(ext$CHROM)

# create list
sub_list <- vector(mode = "list", length = length(ext$CHROM))

for (i in 1:length(ext$CHROM)){
  sub <- res[ext$beg[i] <= res$POS & ext$end[i] >= res$POS & ext$CHROM[i] == res$CHROM]
  sub_list[[i]] <- sub
}

sub <- rbindlist(sub_list)
sub<- sub[,c(1:10)]

#check for dupes
#sub[POS == 4270816 & CHROM == "A07"]

# remove dupes
#sub <- unique(sub, by = c("CHROM", "POS"))
# check again 
#sub[POS == 4270816 & CHROM == "A07"]

# all good so now use this to make the pca 
# write to file then use unix to sort through

sub2<- sub[,c(1:2)]
fwrite(sub2, paste0(subset, "_subset_sigzone_snps.tsv"))

ext2<- ext[,c(1:3)]

# see if things are ordered before
keycol <-c("CHROM","beg")
ext2<- setorderv(ext2, keycol)
colnames(ext2) <- c("chrom", "chromStart","chromEnd")

# need header if using vcftools
fwrite(ext2, sep = "\t", paste0(subset, "_subset_sigzone_snps_vcf.bed"))

## after use script_vcftools_extract
## then run the script_plink_pca_sig_sub

# use this for tabix
fwrite(ext2, sep = "\t", col.names = FALSE, paste0(subset, "_subset_sigzone_snps.bed"))
