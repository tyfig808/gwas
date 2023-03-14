### gapit and then mash
### first convert vcf to hapmap, in bash
PATH=$PATH:/analysis_tools/tassel5.0_standalone
./gwas/tassel-5-standalone/run_pipeline.pl -Xms8G -Xmx8G -fork1 -vcf final_filter.vcf.gz \
    -export final_filter -exportType Hapmap

### in R now
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

### load in phenotypic data
setwd("/net/cephfs/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")
library(tidyverse)

d <- read.csv("phenotypic_data_thomas.csv", stringsAsFactors = T)

### do some cleaning, for full data
d<- d %>% 
  dplyr::rename(
    ID = PlantlabelDNA
  ) %>%
  drop_na(c( 'ID')) 

d<- d %>% 
  subset(d$ID != "RA_TL_H_B_27") #%>%
  #mutate(ID = recode(ID, RA_TL_H_B_27_c = "RA_TL_H_B_27"))

### create subset of pheno traits for gwas, with full traits
p <- d %>% 
  select(c("ID",  "Heightd34", "NbFlowers",  "Indole", "Methyl_anthranilate.1"))

### or with defined subset, lets try with T_H_B because there rep diffs
phen=read.table("pheno_T_H_B.txt", h=T,na.strings=".")
phen <- phen[-6]

gen=read.table("subsetT_H_B_geno.bim", h=T) # a table with first column the chromosomes, and the second column the Position of each marker0 The order should be the same than .map file, seems like we just want col 1 and 4
gen <- gen %>%
        select(c(1,4)) 

colnames(gen)=c("chrom", "pos")

gen <- gen %>% 
        mutate(ID = paste(gen$chrom, gen$pos, sep = '_'))

subsetT_H_B_geno.bed
git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git