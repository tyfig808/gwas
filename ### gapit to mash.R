### gapit to mash
# see vignete here
#https://github.com/Alice-MacQueen/gapit2mashr/blob/master/README.md

### convert gapit output 
### install gapit to mash
#install.packages("devtools")

library(devtools)
#devtools::install_github("Alice-MacQueen/gapit2mashr")
library(gapit2mashr)
library(tidyverse)

setwd("/net/cephfs/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

### redo gapit phenotype function, for use of cmlm
gapit_phenotypes_in_folder <- function(path = ".", model = "CMLM"){
  result_files <- list.files(path = path, pattern = "*GWAS_Results.CMLM*")
  if(model == "CMLM"){
    pheblink <- purrr::partial(stringr::str_sub, start = 37, end = -5)
    # Eventually this needs to accomodate models other than "CMLM", which means
    # modifying where this function ends based on the model used.
    gapit_phenotypes <- purrr::map(result_files, pheblink) %>%
      unlist()
    } else stop("Redo function, use normal for CMLM") 
  return(gapit_phenotypes)
}

### load the phenotypes and check 
phenotypes_vector <- gapit_phenotypes_in_folder(path = "./", model = "CMLM") 
phenotypes_vector

### determine bin of snps to use
### should be ((#SNPS) / (#pheno^2)) = 1,000,000
### if 35 pheno traits, then about 816 snps per bin
numSNPs <- 1000000 / length(phenotypes_vector)^2

path = "./"
## run mashr on file directory
load_GAPIT_GWAS_ty <- function(path = ".", phenotype, model = "CMLM"){
  out <- list()
  out$Pred <- suppressWarnings(readr::read_csv(file.path(path,
                                                         paste0("GAPIT.Association.Pred.", model,
                                                                ".", phenotype,
                                                                ".csv")),
                                               col_names = TRUE,
                                               col_types = "ciiinnnnn"))
  out$Results <- readr::read_csv(file.path(path, paste0("GAPIT.Association.GWAS_Results.", model, ".",
                                                        phenotype,
                                                        ".csv")),
                                 col_names = TRUE, col_types = "ccnnninn") %>%
    dplyr::rename(Chromosome = Chr, Position = Pos, `FDR_Adjusted_P-values` = 'H&B.P.Value') %>%
    dplyr::arrange(.data$Chromosome, .data$Position)

  out$Effects <- readr::read_csv(file.path(path, paste0("GAPIT.Association.GWAS_StdErr.", model, ".",
                                                        phenotype,
                                                      ".csv")),
                                 col_names = TRUE, col_types = "ccninnn") %>%
    dplyr::arrange(.data$Chromosome, .data$Position)
  #out$ROC <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".",
  #                                                  phenotype, ".ROC.csv")),
  #                           col_names = c("Power", "QTN=0.3", "QTN=0.2",
  #                                         "QTN=0.1", "QTN=0.05",
  #                                         "QTN=0.02", "QTN=0.01", "QTN=0"),
  #                           col_types = "nnnnnnnn")
  #Log1c <- readr::read_delim(file.path(path, paste0("GAPIT.", model, ".",
  #                                                  phenotype, ".Log.csv")),
  #                           delim = ";", col_names = "Model_Value",
  #                           col_types = "c")
  #out$Log <- tidyr::separate(Log1c, .data$Model_Value, c("Model", "Value"),
  #                           sep = ",", extra = "merge")
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

gap_to_mash <- function(path = ".", phenotypes = NA, numSNPs = 1000,
                        model = "Blink", S_hat = c("Hedges' G", "ones"),
                        saveoutput = FALSE){
  match.arg(S_hat, c("Hedges' G", "ones"))
  if(is.na(phenotypes)){
    phe_col <- gapit_phenotypes_in_folder(path = path)
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

  df1 <- load_GAPIT_GWAS_ty(path = path, phenotype = phe_col[1])
  big_effects_df <- gapit_top_effects_FDRpvalue(df = df1$Results,
                                                phenotype = phe_col[1],
                                                numSNPs = numSNPs)

  for(i in seq_along(phe_col)[-1]){
    df1 <- load_GAPIT_GWAS_ty(path = path, phenotype = phe_col[i])
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
                                             paste0("effects_", numSNPs,
                                                    "SNPs_PartOneOutput.rds")))
  }

  message(paste0("Part One: data frame of SNPs to keep complete."))
  message(paste0("Starting Part Two: Creating strong and random dataframes of
                 B_hat and S_hat values for use in mashr."))

  df1 <- load_GAPIT_GWAS_ty(path = path, phenotype = phe_col[1])

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

      df1 <- load_GAPIT_GWAS_ty(path = path, phenotype = phe_col[i])

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
gap_to_mash(path = "./", model = "CMLM", numSNPs = numSNPs, S_hat = "Hedges' G", saveoutput = TRUE)

### here is the old function just in case
#gapit2mashr(path = "./", model = "CMLM", numSNPs = numSNPs, S_hat = "Hedges' G", saveoutput = TRUE)

GAPIT.Association.Pred.CMLM.NbFlowers.csv
GAPIT.Association.GWAS_Results.CMLM.Indole.csv