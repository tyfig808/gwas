
## Detecting genomic outliers in a sequence of correlated scores, from Fariello et al 2017.
## Detecting minor and major QTL in GWAS using the local score approach, from Bonhomme et al. 2019.


### Import libraries.
library(ggplot2)
library(data.table)
library(RColorBrewer)

#Import functions.
source('scorelocalfunctions.R')
mydata = fread("example_data_file.txt",h=T)
setkey(mydata, chr)
mydata$pval[mydata$pval==0]=1e-16#min(mydata$pval[-which(mydata$pval ==0)])
Nchr=length(mydata[,unique(chr)])

### Computation of absolute position in the genome. 
#This is useful for doing genomewide plots. 
chrInfo=mydata[,.(L=.N,cor=autocor(pval)),chr]
setkey(chrInfo,chr)
data.table(chr=mydata[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr])))

### Choice of $\xi$ (1,2,3 or 4)

mean(mydata$score)

# The score mean must be negative, ksi must be chosen between mean and max of -log10(p-value)

xi=2
mydata[,score:= -log10(pval)-xi]
mydata[,lindley:=lindley(score),chr]

# verify mean is negative and max is positive
mean(-log10(mydata$pval));max(-log10(mydata$pval))

#hist(mydata$score)
max(mydata$lindley)

# Compute significance threshold for each chromosome
## Uniform distribution of p-values
chrInfo[,th:=thresUnif(L, cor, xi),chr]
(mydata=mydata[chrInfo])

(sigZones=mydata[chrInfo,sig_sl(lindley, pos, unique(th)),chr])


write.table(cbind(mydata$chr,mydata$pos,mydata$pval,mydata$cor,mydata$lindley,mydata$th),file="LocScore_example_data_file_ksi_2.txt",quote=F, sep="\t",row.names=F,col.names=c("chr","pos","pval","cor","lind","th"))
