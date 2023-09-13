library(ggplot2)
library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
all_cov_file <- args[1] # Path to covariance file. Example covariance file see Covariance_example.covar
pheno_name <- args[2] # Name of phenotype
pheno_file <- args[3] # Path to phenotype file for this trait. Example phenotype file see Phenotype_example.pheno
ldpred_filenames <- args[4] # File containing the list of directories of LD Pred results. Each directory should be separated by a line. Methods of results should be in the order of GMMa, GMMb, GMMc, GMMd, MTAG, GWAS, and 4 bins.
pheno_set_name <- args[5] # Phenotype set, named by the initials of each trait in the set. (such as atw, p_ma, p_am, p_ah)
outpath <- args[6]

message(sprintf("PHENO: %s", pheno_name))

all_cov <- read.table(all_cov_file, header = T) %>% select(-c(FID))

pheno <- read.delim(phenofile, header = F)
ldpred_files <- readLines(ldpred_filenames)

###############################################

gmma_lines <- readLines(ldpred_files[1])
gmma_lines <- gsub(", ", ",", gmma_lines)
gmma <- read.table(text=gmma_lines, sep=",", header=T) %>% mutate(PRS_gmma = PRS)

gmmb_lines <- readLines(ldpred_files[2])
gmmb_lines <- gsub(", ", ",", gmmb_lines)
gmmb <- read.table(text=gmmb_lines, sep=",", header=T) %>% mutate(PRS_gmmb = PRS)

gmmc_lines <- readLines(ldpred_files[3])
gmmc_lines <- gsub(", ", ",", gmmc_lines)
gmmc <- read.table(text=gmmc_lines, sep=",", header=T) %>% mutate(PRS_gmmc = PRS)

gmmd_lines <- readLines(ldpred_files[4])
gmmd_lines <- gsub(", ", ",", gmmd_lines)
gmmd <- read.table(text=gmmd_lines, sep=",", header=T) %>% mutate(PRS_gmmd = PRS)

mtag_lines <- readLines(ldpred_files[5])
mtag_lines <- gsub(", ", ",", mtag_lines)
mtag <- read.table(text=mtag_lines, sep=",", header=T) %>% mutate(PRS_mtag = PRS)

gwas_lines <- readLines(ldpred_files[6])
gwas_lines <- gsub(", ", ",", gwas_lines)
gwas <- read.table(text=gwas_lines, sep=",", header=T) %>% mutate(PRS_gwas = PRS)

4bins_lines <- readLines(ldpred_files[7])
4bins_lines <- gsub(", ", ",", 4bins_lines)
4bins <- read.table(text=4bins_lines, sep=",", header=T) %>% mutate(PRS_4bins = PRS)


colnames(pheno) <- c("FID", "IID", pheno_name)

##############################################
merge <- pheno %>% inner_join(all_cov, by = "IID") %>% select(-c(FID))
cov <- merge %>% select(-c(IID))
gmma <- merge %>% inner_join(gmma, by = "IID")  %>% subset(select = -c(IID,true_phens))
gmmb <- merge %>% inner_join(gmmb, by = "IID")  %>% subset(select = -c(IID,true_phens))
gmmc <- merge %>% inner_join(gmmc, by = "IID")  %>% subset(select = -c(IID,true_phens))
gmmd <- merge %>% inner_join(gmmd, by = "IID")  %>% subset(select = -c(IID,true_phens))
mtag <- merge %>% inner_join(mtag, by = "IID")  %>% subset(select = -c(IID,true_phens))
gwas <- merge %>% inner_join(gwas, by = "IID")  %>% subset(select = -c(IID,true_phens))
4bins <- merge %>% inner_join(4bins, by = "IID")  %>% subset(select = -c(IID,true_phens))
8bins <- merge %>% inner_join(8bins, by = "IID")  %>% subset(select = -c(IID,true_phens))


cov_r2 <- summary(lm(percent ~. ,cov))$r.squared
cov_adj_r2 <- summary(lm(percent ~. ,cov))$adj.r.squared

gmma_r2 <- summary(lm(percent ~. ,gmma))$r.squared
gmma_adj_r2 <- summary(lm(percent ~. ,gmma))$adj.r.squared
gmma_inc_r2 <- gmma_r2-cov_r2
gmma_inc_adj_r2 <- gmma_adj_r2-cov_adj_r2

gmmb_r2 <- summary(lm(percent ~. ,gmmb))$r.squared
gmmb_adj_r2 <- summary(lm(percent ~. ,gmmb))$adj.r.squared
gmmb_inc_r2 <- gmmb_r2-cov_r2
gmmb_inc_adj_r2 <- gmmb_adj_r2-cov_adj_r2


gmmc_r2 <- summary(lm(percent ~. ,gmmc))$r.squared
gmmc_adj_r2 <- summary(lm(percent ~. ,gmmc))$adj.r.squared
gmmc_inc_r2 <- gmmc_r2-cov_r2
gmmc_inc_adj_r2 <- gmmc_adj_r2-cov_adj_r2


gmmd_r2 <- summary(lm(percent ~. ,gmmd))$r.squared
gmmd_adj_r2 <- summary(lm(percent ~. ,gmmd))$adj.r.squared
gmmd_inc_r2 <- gmmd_r2-cov_r2
gmmd_inc_adj_r2 <- gmmd_adj_r2-cov_adj_r2


mtag_r2 <- summary(lm(percent ~. ,mtag))$r.squared
mtag_adj_r2 <- summary(lm(percent ~. ,mtag))$adj.r.squared
mtag_inc_r2 <- mtag_r2-cov_r2
mtag_inc_adj_r2 <- mtag_adj_r2-cov_adj_r2


gwas_r2 <- summary(lm(percent ~. ,gwas))$r.squared
gwas_adj_r2 <- summary(lm(percent ~. ,gwas))$adj.r.squared
gwas_inc_r2 <- gwas_r2-cov_r2
gwas_inc_adj_r2 <- gwas_adj_r2-cov_adj_r2


4bins_r2 <- summary(lm(percent ~. ,4bins))$r.squared
4bins_adj_r2 <- summary(lm(percent ~. ,4bins))$adj.r.squared
4bins_inc_r2 <- 4bins_r2-cov_r2
4bins_inc_adj_r2 <- 4bins_adj_r2-cov_adj_r2


8bins_r2 <- summary(lm(percent ~. ,8bins))$r.squared
8bins_adj_r2 <- summary(lm(percent ~. ,8bins))$adj.r.squared
8bins_inc_r2 <- 8bins_r2-cov_r2
8bins_inc_adj_r2 <- 8bins_adj_r2-cov_adj_r2




###############################################
df <- data.frame("Trait" = rep(pheno_name,8), "Method" = c("Cov", "GMMa", "GMMb", "GMMc", "GMMd", "MTAG", "GWAS", "4bins",), "R2" = c(cov_r2, gmma_r2-cov_r2, gmmb_r2-cov_r2, gmmc_r2-cov_r2, gmmd_r2-cov_r2, mtag_r2-cov_r2, gwas_r2-cov_r2, 4bins_r2-cov_r2, 8bins_r2-cov_r2), "Adj_R2" = c(cov_adj_r2, gmma_adj_r2-cov_adj_r2, gmmb_adj_r2-cov_adj_r2, gmmc_adj_r2-cov_adj_r2, gmmd_adj_r2-cov_adj_r2, mtag_adj_r2-cov_adj_r2, gwas_adj_r2-cov_adj_r2, 4bins_adj_r2-cov_adj_r2, 8bins_adj_r2-cov_adj_r2))

df

saveRDS(df, file = file.path(outpath, sprintf("%s_df_%s.rds", pheno_name, pheno_set_name)))
