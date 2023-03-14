### no run mash after running gapit, and gapit to mash
### see vignete here
#https://github.com/Alice-MacQueen/CDBNgenomics/blob/master/README.md

### check to see if need additional packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

### if want to install packages for annotation 
BiocManager::install(c("multtest", "GenomicFeatures", "GenomicRanges", 
                     "IRanges", "VariantAnnotation", "AnnotationDbi"))

library(devtools)
#devtools::install_github("Alice-MacQueen/CDBNgenomics")

# install mashr
#devtools::install_github("stephenslab/mashr@v0.2-11")

# for plotting
#devtools::install_github("lcolladotor/dots")


library(CDBNgenomics)
library(mashr)
library(ashr)
library(gapit2mashr)
library(tidyverse)

setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

### run mash, if just after runnng the gapit to mash
mash_output <- mash_standard_run(path = ".", list_input = mash_input, 
                                  numSNPs = numSNPs, saveoutput = TRUE)


# Or, if you are doing this in a new session and don't have mash_input 
#     in your workspace, you just need to enter the number of SNPs and this 
#     function will find a previously saved rds file with that numSNPs for you:
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
phenotypes_vector <- gapit_phenotypes_in_folder(path = "./", model = "CMLM") 
numSNPs <- 1000000 / length(phenotypes_vector)^2

mash_output <- mash_standard_run(path = ".", numSNPs = numSNPs, saveoutput = TRUE)
# if get error check this
# so it put this out 
#effects_62500SNPs_PartOneOutput.rds

### we want this
#Part-One-Output_Top-Effects-62500-SNPs.rds

### visualize mash
### manhattan with pleio

mashhattan <- mash_plot_manhattan_by_condition(m = mash_output)
#> Joining, by = "value"
mashhattan$ggmanobject

pairwise_plot <- mash_plot_pairwise_sharing(effectRDS = "Pairwise_sharing_Strong_Effects_62500SNPs.rds",
                                            reorder = TRUE) 
#> Loading required namespace: dots
#> Scale for 'colour' is already present. Adding another scale for 'colour',
#> which will replace the existing scale.
pairwise_plot$gg_corr

### effects plot
effects <- mash_plot_effects(m = mash_output, n = 1)
effects$ggobject +
    scale_x_discrete(labels = str_sub(as_vector(effects$effect_df$value), 
                                      start = 6)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
#> Warning: Unknown or uninitialised column: `value`.