## mash from emmax, try to see pleiotropy
# check nice vignette here
# https://github.com/stephenslab/mashr/blob/master/vignettes/intro_mash.Rmd

# load packages and wd
library(ashr)
library(mashr)
library(vroom)
library(CDBNgenomics) # can use this to access her plots


setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

# read in data
#Bhat <- vroom("B_hat_random_df_62500topSNPs")

# if loading in from gapit to mash
B_hat_random_df_62500topSNPs <- readRDS("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas/B_hat_random_df_62500topSNPs.rds")
S_hat_random_df_62500topSNPs <- readRDS("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas/S_hat_random_df_62500topSNPs.rds")

# load in bhat strong, for some reason, these are only with ones from chrom 1, need to check the gapit to mash code
B_hat_strong_df_62500topSNPs <- readRDS("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas/B_hat_strong_df_62500topSNPs.rds")
S_hat_strong_df_62500topSNPs <- readRDS("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas/S_hat_strong_df_62500topSNPs.rds")
# tells me it needs something else for type list
bhat <- as.matrix(B_hat_random_df_62500topSNPs)
shat <- as.matrix(S_hat_random_df_62500topSNPs)

bhat_strong<- as.matrix(B_hat_strong_df_62500topSNPs)
shat_strong <- as.matrix(S_hat_strong_df_62500topSNPs)
# need bhat and shat (effect and std error)
data = mash_set_data(bhat, shat)
data_strong = mash_set_data(bhat_strong, shat_strong)

# set up covariance matrics, pc must be more than n_conditions(data)
U.pca = cov_pca(data_strong,3) # create pca from strong dataset, data driven covariate
U.ed = cov_ed(data_strong, U.pca) # this can take a while, less than 50 min for 62500 snps

# or use this canonical matrix
U.c = cov_canonical(data)  
print(names(U.c))

# run the model, use mash on both covarinace structures
m = mash(data, Ulist = c(U.ed,U.c)) 

# here we run the model than update it to use it on the strong effects
m = mash(data, Ulist = c(U.ed,U.c), outputlevel = 1)
m2 = mash(data_strong, g=get_fitted_g(m), fixg=TRUE)  

# get posterior means, can convert these into pvalues, 
# lea method is to take rank of posterior mean and then divide by number of snps to get arbitrary pvalue
head(get_pm(m))

# get sig results, and see how many
head(get_significant_results(m))
print(length(get_significant_results(m)))

# to get results in subset of condtion, not sure what cond = 1 refers to, hopefully the first trait
print(head(get_significant_results(m, conditions=1)))

mashhattan <- mash_plot_manhattan_by_condition(m = m)
#> Joining, by = "value"
mashhattan$ggmanobject


mashhattan2 <- mash_plot_manhattan_by_condition(m = m2)
#> Joining, by = "value"
mashhattan2$ggmanobject
# get pairwise sharing, need to check how to use function where only look at negative matches to get antonistic pleio
print(get_pairwise_sharing(m)) 
print(get_pairwise_sharing(m, FUN=abs)) # sharing by magnitude when sign is ignored

# to see where most of model variation is
print(get_estimated_pi(m))
barplot(get_estimated_pi(m),las = 2)

# if ever get error that plot margins are too large
par(mar=c(9,9,9,9))

# get number of sig conditions per snp
x <- get_n_significant_conditions(m)
x <- x %>% 
	separate(x, c("snp", "# cond"), sep = "\t") 

# effects plots
effects <- mash_plot_effects(m = m, n = 1)
effects$ggobject +
    scale_x_discrete(labels = str_sub(as_vector(effects$effect_df$value), 
                                      start = 6)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))