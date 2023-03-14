### emmax, lindley score to mash
## so here we want to conver the lindley score as the adjust pvalue 
## to a mash data frame, emmax provides the effect (hhat)
## and also the standard error of the effect (shat)
# see vignete here
#https://github.com/Alice-MacQueen/gapit2mashr/blob/master/README.md

### install gapit to mash
#install.packages("devtools")

library(devtools)
#devtools::install_github("Alice-MacQueen/gapit2mashr")
library(gapit2mashr)
library(tidyverse)
library(vroom)
## use this for rstudio, need another for the cluster
setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

### here enter the data subset for usage
subset <- "T_H_B"

### now read in phenotype.txt needed for emmeans
### then get the number of pheno to determine number of snps
phen=read.table(paste0("pheno_", subset,".txt"), h=T,na.strings=".")
phen<-phen[1:(length(phen)-1)]

phenotypes_vector=colnames(phen)[2:ncol(phen)]
numSNPs <- 1000000 / length(phenotypes_vector)^2

path = "."


### so need to load the emmax.ps, snp, beta, std error, and pvalue
## check the files of gapit to see, but also need maf somewhere, could grab from the freq file

#phenotype = phenotypes_vector[1]

### for the results file we need freq, emmax, and the lindley file
freq <- vroom(file.path(path, paste0("subset", subset, "_drift.frq")), delim = "\t", 
      col_names = c("Chromosome", "Position", "nalleles", "nchr", "a1", "a2"), skip = 1) %>%
      dplyr::mutate(SNP = paste(Chromosome, Position, sep = '_')) %>%
      separate(a2, c("a2_base", "maf"), sep = ":")

freq$maf <- as.numeric(freq$maf)
# made a mistake with the mrk other files so remove from freq and effect
freq = freq[-1,]

### load in emmax effect file
emmax <- vroom(file.path(path,paste0(subset, "_", phenotype, ".emmax.ps")),
           col_names = c("SNP", "effect", "std Error", "pvalue")) %>%
    dplyr::select(-1) 
emmax = emmax[-1,] # redo for the others so dont have to drop

### load in lindley file
lindley <- vroom(file.path(path, paste0("Local_Score_", phenotype, 
                                                        "_xi2.txt"))) %>%
    ## here we use the lindley instead of the fdr value but we just call it for easy code
    dplyr::rename(Chromosome = chrom, Position = POS, `FDR_Adjusted_P-values` = 'lindley') %>%
    dplyr::arrange(.data$Chromosome, .data$Position) %>% 
    dplyr::mutate(SNP = paste(Chromosome, Position, sep = '_'))

### combine into results file
Results <- bind_cols(freq$SNP, freq$Chromosome, freq$Position, freq$maf, 
                     emmax$pvalue, lindley$`FDR_Adjusted_P-values`, emmax$effect) %>% 
        dplyr::rename(
        SNP = ...1,
        Chromosome = ...2,
        Position =  ...3,
        maf =  ...4,
        pvalue =  ...5,
        `FDR_Adjusted_P-values` =  ...6,
        effect =  ...7) %>%
        dplyr::arrange(.data$Chromosome, .data$Position) 

x <- get(paste0("phen"))[[phenotype]] 
Results$nobs <- sum(!is.na(x))-1

### combine into effects file
Effects <- bind_cols(freq$SNP, freq$Chromosome, freq$Position,  
                     emmax$'std Error', emmax$effect) %>% 
        dplyr::rename(
        SNP = ...1,
        Chromosome = ...2,
        Position =  ...3,
        'std Error' =  ...4,
        effect =  ...5) %>%
        dplyr::arrange(.data$Chromosome, .data$Position) 

# need degrees of freedom, 
# maybe we assume for emmax that it is the number of obs - 1 
x <- get(paste0("phen"))[[phenotype]] 
Effects$df <- sum(!is.na(x))-2

### new function with emmax, beacuse the trait is in the middle has trouble grabbing, try later 
gapit_phenotypes_in_folder <- function(path = ".", subset = "subset", model = "emmax"){
  p<-paste0('"*', subset, "_", phenotype, '.emmax.ps*"')
  result_files <- list.files(path = path, pattern = noquote(p))
  if(model == "emmax"){
    pheblink <- purrr::partial(stringr::str_sub, start = 7, end = -5)
    # Eventually this needs to accomodate models other than "CMLM", which means
    # modifying where this function ends based on the model used.
    gapit_phenotypes <- purrr::map(result_files, pheblink) %>%
      unlist()
    } else stop("Redo function, use normal for CMLM") 
  return(gapit_phenotypes)
}

# for testing
#y <- gapit_phenotypes_in_folder(path = ".", subset = subset)
#df1 <- load_emmax_GWAS_ty(path = path, phenotype = "Heightd34", subset = subset)
### try function
load_emmax_GWAS_ty <- function(path = ".", phenotype, subset = "subset", freq = freq){
  out <- list()
  # load emmax
  emmax <- vroom(file.path(path,paste0(subset, "_", phenotype, ".emmax.ps")),
           col_names = c("SNP", "effect", "std Error", "pvalue")) %>%
    dplyr::select(-1) 
  emmax = emmax[-1,] # redo for the others so dont have to drop
  # load lindley
  lindley <- vroom(file.path(path, paste0("Local_Score_", phenotype, 
                                                        "_xi2.txt"))) %>%
    ## here we use the lindley instead of the fdr value but we just call it for easy code
    dplyr::rename(Chromosome = chrom, Position = POS, `FDR_Adjusted_P-values` = 'lindley') %>%
    dplyr::arrange(.data$Chromosome, .data$Position) %>% 
    dplyr::mutate(SNP = paste(Chromosome, Position, sep = '_'))
  # now create outfile for later functions
  out$Results <- bind_cols(freq$SNP, freq$Chromosome, freq$Position, freq$maf, 
                     emmax$pvalue, lindley$`FDR_Adjusted_P-values`, emmax$effect) %>% 
        dplyr::rename(
        SNP = ...1,
        Chromosome = ...2,
        Position =  ...3,
        maf =  ...4,
        pvalue =  ...5,
        `FDR_Adjusted_P-values` =  ...6,
        effect =  ...7) %>%
        dplyr::arrange(.data$Chromosome, .data$Position) 
  # get phenotype trait, then get non na
  x <- get(paste0("phen"))[[phenotype]] 
  out$Results$nobs <- sum(!is.na(x))-1
  # combine into effects file
  out$Effects <- bind_cols(freq$SNP, freq$Chromosome, freq$Position,  
                     emmax$'std Error', emmax$effect) %>% 
        dplyr::rename(
        SNP = ...1,
        Chromosome = ...2,
        Position =  ...3,
        'std Error' =  ...4,
        effect =  ...5) %>%
        dplyr::arrange(.data$Chromosome, .data$Position) 

  # maybe we assume for emmax that it is the number of obs - 1 
  out$Effects$df <- sum(!is.na(x))-2

  return(out)
}

gapit_top_effects_FDRpvalue <- function(df, phenotype, numSNPs){
  dfA <- df %>%
    dplyr::arrange(.data$`FDR_Adjusted_P-values`) %>%
    dplyr::select(-tidyselect::starts_with("Rsquare"))
  df2 <- dfA[1:numSNPs,]
  names(df2)[4] <- paste0(phenotype, "_pvalue")
  names(df2)[5] <- paste0(phenotype, "_maf")
  names(df2)[6] <- paste0(phenotype, "_nobs")
  names(df2)[7] <- paste0(phenotype, "_FDR_adj_pvalue")
  names(df2)[8] <- paste0(phenotype, "_effect")
  return(df2)
}

s_hat_hedges_g <- function(df, phenotype){
  standardization <- max(abs(df$effect), na.rm = TRUE)
  df3 <- df %>%
    dplyr::mutate(Stand_effect = .data$effect / standardization,
                  Obs = .data$maf * .data$nobs,
                  Obs2 = (1-.data$maf) * .data$nobs,
                  d = ifelse(abs(.data$Stand_effect) < 0.98,
                             (2 * .data$Stand_effect) /
                               sqrt(1 - .data$Stand_effect^2),
                             4),
                  d_unbiased = (1 - (3 / (4 * (.data$nobs -2) -1))) * .data$d,
                  sigma2_d = ((.data$Obs + .data$Obs2) /
                                (.data$Obs * .data$Obs2)) +
                    (.data$d_unbiased^2 / (2*(.data$Obs + .data$Obs2))),
                  stderr_d = sqrt(.data$sigma2_d)) %>%
    dplyr::mutate(Stand_effect = ifelse(is.na(.data$Stand_effect) |
                                          is.infinite(.data$Stand_effect),
                                        0,
                                        .data$Stand_effect),
                  stderr_d = ifelse(is.na(.data$stderr_d) |
                                      is.infinite(.data$stderr_d),
                                    10,
                                    .data$stderr_d))  %>%
    dplyr::select(.data$SNP, .data$Stand_effect, .data$stderr_d)
  names(df3)[2] <- paste0("Bhat_", phenotype)
  names(df3)[3] <- paste0("Shat_", phenotype)
  return(df3)
}

s_hat_gapit <- function(df, phenotype){
  standardization <- max(abs(df$effect), na.rm = TRUE)

  df3 <- df %>%
    dplyr::mutate(stderr_d = .data$`std Error` / standardization,
                  Stand_effect = .data$effect / standardization) %>%
    dplyr::select(.data$SNP, .data$Stand_effect, .data$stderr_d) %>%
    dplyr::mutate(stderr_d = ifelse(is.na(.data$stderr_d),
                                    10,
                                    .data$stderr_d),
                  Stand_effect = ifelse(is.na(.data$Stand_effect),
                                        0,
                                        .data$Stand_effect))
  names(df3)[2] <- paste0("Bhat_", phenotype)
  names(df3)[3] <- paste0("Shat_", phenotype)
  return(df3)
}

# for testing
#gapit_top_effects_FDRpvalue(x, phenotype = NA)
#phe_col <- gapit_phenotypes_in_folder(path = path)


emmax_to_mash <- function(path = ".", phenotypes = NA, numSNPs = 1000,
                        subset = "subset", S_hat = c("Hedges' G", "ones"),
                        saveoutput = FALSE){
  match.arg(S_hat, c("Hedges' G", "ones"))
  if(is.na(phenotypes)){
    phen=read.table(paste0("pheno_", subset,".txt"), h=T,na.strings=".")
    # drop the last column
    phen<-phen[1:(length(phen)-1)]
    phe_col <- colnames(phen)[2:ncol(phen)]
  } else {
    phe_col <- phenotypes
  }
  if(is.null(phe_col)){
    stop("Can't find any GAPIT Results files in this path.")
  }
  if(is.na(phe_col[1])){
    stop("Can't find any GAPIT Results files in this path.")
  }
  numSNPs <- as.numeric(numSNPs)
  message(paste0("Starting part one: Making a data frame of all SNPs that are",
                 " in the top ", numSNPs, " SNPs
                 by FDR adjusted p-values for at least one phenotype."))
  # load in freq so dont need to load in everytime
  freq <- vroom(file.path(path, paste0("subset", subset, "_drift.frq")), delim = "\t", 
      col_names = c("Chromosome", "Position", "nalleles", "nchr", "a1", "a2"), skip = 1) %>%
      dplyr::mutate(SNP = paste(Chromosome, Position, sep = '_')) %>%
      separate(a2, c("a2_base", "maf"), sep = ":")
  freq$maf <- as.numeric(freq$maf)
  # made a mistake with the mrk other files so remove from freq and effect
  freq = freq[-1,]

  df1 <- load_emmax_GWAS_ty(path = path, phenotype = phe_col[1], subset = subset, freq = freq)
  big_effects_df <- gapit_top_effects_FDRpvalue(df = df1$Results,
                                                phenotype = phe_col[1],
                                                numSNPs = numSNPs)

  for(i in seq_along(phe_col)[-1]){
    df1 <- load_emmax_GWAS_ty(path = path, phenotype = phe_col[i], subset = subset, freq = freq)
    df2 <- gapit_top_effects_FDRpvalue(df = df1$Results,
                                       phenotype = phe_col[i],
                                       numSNPs = numSNPs)
    big_effects_df <- df2 %>%
      dplyr::full_join(big_effects_df, by = c("SNP", "Chromosome",
                                              "Position"))
  }
  big_effects_df <- big_effects_df %>%
    dplyr::arrange(.data$Chromosome, .data$Position)

  if(saveoutput == TRUE){
    saveRDS(big_effects_df, file = file.path(path,
                                             paste0("Part-One-Output_Top-Effects-", numSNPs,
                                                    "-SNPs.rds")))
  }

  message(paste0("Part One: data frame of SNPs to keep complete."))
  message(paste0("Starting Part Two: Creating strong and random dataframes of
                 B_hat and S_hat values for use in mashr."))

  df1 <- load_emmax_GWAS_ty(path = path, phenotype = phe_col[1], subset = subset, freq = freq)

  if(S_hat == "Hedges' G"){ # fix this: need support for "ones"
    if(sum(is.na(df1$Effects$`std Error`)) > length(df1$Effects$`std Error`)*.05){
      # If there are too many NA's for standard errors, derive new standard errors
      # using Hedges' G (which requires the MAF).
      df3 <- s_hat_hedges_g(df = df1$Results, phenotype = phe_col[1])
      message(paste0("Hedge's G standard errors were used for the phenotype '",
                     phe_col[1], "'."))
    } else {
      # or if not many of the standard errors are NA's, just use them for Shats.
      df3 <- s_hat_gapit(df = df1$Effects, phenotype = phe_col[1])
      message(paste0("GAPIT's standard errors were used for the phenotype '",
                     phe_col[1], "'."))
    }

    # Start making data frames of strong and random B_hat and S_hat.
    bhat_df <- big_effects_df %>%
      dplyr::select(.data$SNP) %>%
      dplyr::left_join(df3, by = "SNP") %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_df <- big_effects_df %>%
      dplyr::select(.data$SNP) %>%
      dplyr::left_join(df3, by = "SNP") %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))

    set.seed(1234) # Makes the random data frames reproducible.
    random_sample <- sample(1:nrow(df3), nrow(big_effects_df)) %>% sort()
    bhat_random <- df3[random_sample,] %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_random <- df3[random_sample,] %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))

    for(i in seq_along(phe_col)[-1]){

      df1 <- load_emmax_GWAS_ty(path = path, phenotype = phe_col[i], subset = subset, freq = freq)

      if(sum(is.na(df1$Effects$`std Error`)) > length(df1$Effects$`std Error`)*.05){
        # If there are too many NA's for standard errors, derive new standard
        # errors using Hedge's G (which requires the MAF).
        df3 <- s_hat_hedges_g(df = df1$Results, phenotype = phe_col[i])
        message(paste0("Hedge's G standard errors were used for the phenotype '",
                       phe_col[i], "'."))
      } else {
        # if not many of the standard errors are NA's, just use them for Shats.
        df3 <- s_hat_gapit(df = df1$Effects, phenotype = phe_col[i])
        message(paste0("GAPIT's standard errors were used for the phenotype '",
                       phe_col[i], "'."))
      }

      bhat_df <- bhat_df %>%
        dplyr::left_join(df3, by = "SNP") %>%
        dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
      shat_df <- shat_df %>%
        dplyr::left_join(df3, by = "SNP") %>%
        dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
      bhat_random <- bhat_random %>%
        dplyr::left_join(df3[random_sample,], by = "SNP") %>%
        dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
      shat_random <- shat_random %>%
        dplyr::left_join(df3[random_sample,], by = "SNP") %>%
        dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
    }
  }

  B_hat_random <- data.frame(bhat_random, row.names = "SNP")
  S_hat_random <- data.frame(shat_random, row.names = "SNP")
  B_hat_strong <- data.frame(bhat_df, row.names = "SNP")
  S_hat_strong <- data.frame(shat_df, row.names = "SNP")

  if(saveoutput == TRUE){
    saveRDS(B_hat_strong, file = file.path(path, paste0("B_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_strong, file = file.path(path, paste0("S_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(B_hat_random, file = file.path(path, paste0("B_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_random, file = file.path(path, paste0("S_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
  }
  return(list(SNP_df = big_effects_df, B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong, B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}

### try your function
emmax_to_mash(path = path, numSNPs = numSNPs, subset = subset, S_hat = "Hedges' G", saveoutput = TRUE)

### here is the old function just in case
#gapit2mashr(path = "./", model = "CMLM", numSNPs = numSNPs, S_hat = "Hedges' G", saveoutput = TRUE)
