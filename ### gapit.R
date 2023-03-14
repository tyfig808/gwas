### gapit gwas
setwd("/net/cephfs/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

#if (!require(devtools)) install.packages('devtools')
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)

if (!require(vroom)) install.packages('vroom')
if (!require(tidyverse)) install.packages('tidyverse')

library(tidyverse)
library(GAPIT3)
library(vroom)

### for gapit we need phen, and either hapmap or ped/map
### need to call as my_ or else graphing functions dont work down the line
### load phen
myY = vroom("pheno_T_H_B.txt")
myY <- myY[-6]
myY = as.data.frame(myY)
### in bash, vcf to hapmap, takes a while but seems to be easier on r to handle
#./tassel-5-standalone/run_pipeline.pl -Xms8G -Xmx8G -fork1 -vcf subsetT_H_B.recode.vcf \
#    -export subset_T_H_B -exportType Hapmap

### so getting error that x is not numeric, some guy said to do this
myG <- read.delim("./subset_T_H_B.hmp.txt", head = FALSE)

### load in kin file
myKI <- read.table("./T_H_B.aIBS.kinf")
myKI <- cbind(myKI, myY$IID)
myKI <- myKI %>%
        dplyr::rename(ID = `myY$IID`) %>%
        relocate(ID)

## run model
myGAPIT=GAPIT(
Y=myY, #fist column is ID
G=myG,
KI=myKI,
PCA.total=3,
model=c("FarmCPU", "Blink"),
Multiple_analysis=TRUE)
