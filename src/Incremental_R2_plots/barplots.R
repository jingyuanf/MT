library(ggplot2)
library(tidyverse)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
rds_filenames <- args[1] # File containing names of all rds files produced in generate_r2_df.R for phenotypes in a set of analysis
pheno_names <- args[2] # File containing names of all of those phenotypes
color_names <- args[3] # File containing color codes for all phenotypes (number of color codes should match the number of phenotypes)
set_name <- args[4] # The name of the trait set. For example, for armfat_percent, trunkfat_percent and waistc, the name of the trait set could be "atw"
output_path <- args[5] # The output path of the plots.

rdsL <- readLines(rds_filenames)
pheno_nameL <- readLines(pheno_names)
color_nameL <- readLines(color_names)

set <- readRDS(rdsL[1])
set_cov <- filter(set, Method == "Cov")
set$R2 <- set$R2 / set_cov$R2
set$Adj_R2 <- set$Adj_R2 / set_cov$Adj_R2
set <- set %>% filter(Method != "Cov")

for(rds_filename in rdsL[2:length(rdsL)]){
	new_trait <- readRDS(rds_filename)
	new_cov <- filter(new_trait, Method == "Cov")
	new_trait$R2 <- new_trait$R2 / new_cov$R2
	new_trait$Adj_R2 <- new_trait$Adj_R2 / new_cov$Adj_R2
	new_trait <- new_trait %>% filter(Method != "Cov")
	set <- rbind(set, new_trait)
}


cols <- setNames(data.frame(matrix(ncol = 3, nrow = 1)), pheno_nameL)

for (i in 1:length(pheno_nameL)){
	cols[pheno_nameL[i]] <- color_nameL[i]
}

print(cols)


ggplot(set, aes(Method, R2, fill = Trait)) +
  geom_bar(stat = "identity", position = 'dodge') +
  theme_bw()+
  theme(axis.text=element_text(size=20),
	axis.text.x=element_text(angle=60),
        axis.title=element_text(size=20),
	plot.title = element_text(size=25), 
        strip.text = element_text(size=20)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  scale_colour_manual(
    values = cols,
    aesthetics = c("fill")
  ) +
  labs(title="Incremental R2 by Methods",
         x="Trait",
         y="Incremental R2") +
  facet_wrap(~ Trait) +
ggsave(file.path(output_path, sprintf("%s_r2.png", set_name)), width = 40, height = 28, unit = "cm")

ggplot(set, aes(Method, Adj_R2, fill = Trait)) +
  geom_bar(stat = "identity", position = 'dodge') +
  theme_bw()+
  theme(axis.text=element_text(size=20),
	axis.text.x=element_text(angle=60),
        axis.title=element_text(size=30),
	plot.title = element_text(size=25), 
        strip.text = element_text(size=20)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  labs(title="Incremental Adjusted R2 by Methods",
         x="Trait",
         y="Incremental Adjusted R2") +
  scale_colour_manual(
    values = cols,
    aesthetics = c("fill")
  ) +
  facet_wrap(~ Trait) +
ggsave(file.path(output_path, sprintf("%s_adj_r2.png", set_name)), width = 40, height = 28, unit = "cm")

