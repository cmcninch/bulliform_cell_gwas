### Load libraries and set working directory###
##################################################
library(tidyverse)
library(data.table)
library(lubridate)
library(lme4)
library(qvalue)
library(asbio)
library(broom)

setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Flowering_Time_GWAS_Results/")

### Create functions ###
##################################################
blup_single_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Rep) + (1 | SNP_Genotype_Name:Rep), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "Taxa")
  colnames(blups) <- c("Taxa", trait)
  return(blups)
}

### Compute BLUPs for effect size percent estimates ###
##################################################
### Trait enntry-means ###
flowering_times_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv") %>%
  mutate(DTA = as.numeric(mdy(Anthesis_Date) - mdy("5/8/17")), 
         DTS = as.numeric(mdy(Silking_Date) - mdy("5/8/17"))) %>%
  dplyr::select(SNP_Genotype_Name, Row, Rep, Block, DTA, DTS) %>%
  gather("Trait", "Score", 5:6) %>%
  group_by(SNP_Genotype_Name, Row, Rep, Block, Trait) %>%
  summarise(Score = mean(Score, na.rm = TRUE)) %>%
  spread(Trait, Score)

### BLUPs ###
dta_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTA")
dts_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTS")

combined_blups <- full_join(dta_blups, dts_blups, by = "Taxa")

blup_sds <- tibble(Trait = c("DTA", "DTS"),
                   SD = c(sd(dta_blups$DTA, na.rm = TRUE),
                          sd(dts_blups$DTS, na.rm = TRUE)))

### Compile significant SNPs based on Bonferroni threshold ###
##################################################
path <- getwd()
sig_snps <- tibble()
file.names <- dir(path, pattern =".GWAS.csv")
for(i in 1:length(file.names)){
  file <- read_csv(file.names[i], 
                   col_names = TRUE)
  file <- file %>% 
    mutate(Trait = gsub("\\.GWAS\\.csv", "", file.names[i]), 
           `-log10p.value` = as.numeric(-log10(p.value)))

  file$q_value <- qvalue(file$p.value)$qvalues
  
  file <- filter(file, q_value < 0.05)
  
  sig_snps <- bind_rows(sig_snps, file)
}

sig_snps <- dplyr::select(sig_snps, Trait, SNP, `-log10p.value`, q_value, estimate) %>%
  left_join(., blup_sds, by = "Trait") %>%
  mutate(Standardized_Effect_Size = 2*estimate/SD) %>%
  dplyr::select(-estimate, -SD) %>%
  mutate(`-log10p.value` = round(`-log10p.value`, digits = 2), 
         `q_value` = formatC(`q_value`, format = "e", digits = 2), 
         Standardized_Effect_Size = round(Standardized_Effect_Size, digits = 2))

### Add major/minor alleles and minor allele frequencies ###
##################################################
sig_snps <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/snp_general_info.txt") %>%
  left_join(sig_snps, ., by = c("SNP" = "Site Name")) %>%
  mutate(`Major/Minor Allele` = paste(`Major Allele`, "/", `Minor Allele`, sep = "")) %>%
  dplyr::select(Trait, SNP, `-log10p.value`, q_value, `Major/Minor Allele`, `Minor Allele Frequency`, Standardized_Effect_Size)

### Add proportion of genotypic variance explained for each snp; partial R2 ###
##################################################
snps <- fread("../snps/MAF_MissRate_Filtered_WiDiv_SNPs.xmat")

### % of additional variance explained by utilizing the SNP ###
### DTA ###
dta_sig_snps <- as.data.frame(snps)[, c(1, which(colnames(snps) %in% filter(sig_snps, Trait == "DTA")$SNP))] %>%
  left_join(dta_blups, ., by = c("Taxa" = "taxa"))

dta_partial_R2 <- tibble()
snp_list <- colnames(dta_sig_snps)[-c(1:2)]
for(i in 1:length(snp_list)){
  dta_full <- lm(DTA ~ ., data = dta_sig_snps[-1])
  dta_partial <- lm(as.formula(paste("DTA~", paste(snp_list[-i], collapse = "+"))), data = dta_sig_snps[-1])
  temp_partial_R2 <- tibble(SNP = snp_list[i], 
                            Partial_R2 = round(partial.R2(dta_partial, dta_full), digits = 3))
  dta_partial_R2 <- bind_rows(dta_partial_R2, temp_partial_R2)
}
dta_partial_R2 <- mutate(dta_partial_R2, Trait = "DTA")

### DTS ###
dts_sig_snps <- as.data.frame(snps)[, c(1, which(colnames(snps) %in% filter(sig_snps, Trait == "DTS")$SNP))] %>%
  left_join(dts_blups, ., by = c("Taxa" = "taxa"))

dts_partial_R2 <- tibble()
snp_list <- colnames(dts_sig_snps)[-c(1:2)]
for(i in 1:length(snp_list)){
  dts_full <- lm(DTS ~ ., data = dts_sig_snps[-1])
  dts_partial <- lm(as.formula(paste("DTS~", paste(snp_list[-i], collapse = "+"))), data = dts_sig_snps[-1])
  temp_partial_R2 <- tibble(SNP = snp_list[i], 
                            Partial_R2 = round(partial.R2(dts_partial, dts_full), digits = 3))
  dts_partial_R2 <- bind_rows(dts_partial_R2, temp_partial_R2)
}
dts_partial_R2 <- mutate(dts_partial_R2, Trait = "DTS")

partial_R2s <- bind_rows(dta_partial_R2, dts_partial_R2)

### % of variance explained by each SNP when all SNPs in model ###
### DTA ###
dta_aov_partial_r2 <- aov(DTA ~ ., data = dplyr::select(dta_sig_snps, -Taxa)) %>%
  tidy() %>%
  mutate(partial_r2_all = sumsq/sum(sumsq), 
         Trait = "DTA") %>%
  dplyr::select(Trait, term, partial_r2_all) %>%
  dplyr::rename("SNP" = term)

### DTS ###
dts_aov_partial_r2 <- aov(DTS ~ ., data = dplyr::select(dts_sig_snps, -Taxa)) %>%
  tidy() %>%
  mutate(partial_r2_all = sumsq/sum(sumsq), 
         Trait = "DTS") %>%
  dplyr::select(Trait, term, partial_r2_all) %>%
  dplyr::rename("SNP" = term)

aov_partial_r2 <- bind_rows(dta_aov_partial_r2, dts_aov_partial_r2) %>%
  mutate(partial_r2_all = round(partial_r2_all, digits = 3))

### Create final table ###
sig_snps <- left_join(sig_snps, partial_R2s, by = c("Trait", "SNP")) %>%
  left_join(., aov_partial_r2, by = c("Trait", "SNP"))
