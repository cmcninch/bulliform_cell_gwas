### Load libraries and set working directory###
##################################################
library(tidyverse)
library(data.table)
library(lubridate)
library(lme4)
library(qvalue)
library(broom)
library(asbio)

setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Bulliform_GWAS_Results/")

### Create functions ###
##################################################
blup_multiple_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Location) + (1 | SNP_Genotype_Name:Location) + (1 | Rep:Location), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "Taxa") 
  colnames(blups) <- c("Taxa", trait)
  return(blups)
}

### Compute BLUPs for effect size percent estimates ###
##################################################
### Trait enntry-means ###
bulliform_ff_entry_means <- fread("../trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  mutate(Bulliform_FF = Bulliform_Cell_Field_Count/1040 * 245) %>% # height of image is 1040 pixels; 245 pixels = 1mm
  group_by(SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FF = mean(Bulliform_FF, na.rm = TRUE)) %>%
  drop_na()

bulliform_fw_entry_means <- fread("../trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  group_by(Image, SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_Field_Area, na.rm = TRUE)/1388/950*1000) %>% # width of image is 1388 pixels; 950 pixels = 1mm
  group_by(SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_FW, na.rm = TRUE)) %>%
  drop_na()

### BLUPs ###
bulliform_ff_blups <- blup_multiple_enviro_func(df = bulliform_ff_entry_means, trait = "Bulliform_FF")
bulliform_fw_blups <- blup_multiple_enviro_func(df = bulliform_fw_entry_means, trait = "Bulliform_FW")

blup_sds <- tibble(Trait = c("BFF", "BFW"),
                   SD = c(sd(bulliform_ff_blups$Bulliform_FF, na.rm = TRUE),
                          sd(bulliform_fw_blups$Bulliform_FW, na.rm = TRUE)))

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
  mutate(Trait = ifelse(Trait == "Bulliform_FF", "BFF", "BFW")) %>%
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
### BFF ###
bff_sig_snps <- as.data.frame(snps)[, c(1, which(colnames(snps) %in% filter(sig_snps, Trait == "BFF")$SNP))] %>%
  left_join(bulliform_ff_blups, ., by = c("Taxa" = "taxa"))

bff_partial_R2 <- tibble()
snp_list <- colnames(bff_sig_snps)[-c(1:2)]
for(i in 1:length(snp_list)){
  bff_full <- lm(Bulliform_FF ~ ., data = bff_sig_snps[-1])
  bff_partial <- lm(as.formula(paste("Bulliform_FF~", paste(snp_list[-i], collapse = "+"))), data = bff_sig_snps[-1])
  temp_partial_R2 <- tibble(SNP = snp_list[i], 
                            Partial_R2 = round(partial.R2(bff_partial, bff_full), digits = 3))
  bff_partial_R2 <- bind_rows(bff_partial_R2, temp_partial_R2)
}
bff_partial_R2 <- mutate(bff_partial_R2, Trait = "BFF")

### BFW ###
bfw_sig_snps <- as.data.frame(snps)[, c(1, which(colnames(snps) %in% filter(sig_snps, Trait == "BFW")$SNP))] %>%
  left_join(bulliform_fw_blups, ., by = c("Taxa" = "taxa"))

bfw_partial_R2 <- tibble()
snp_list <- colnames(bfw_sig_snps)[-c(1:2)]
for(i in 1:length(snp_list)){
  bfw_full <- lm(Bulliform_FW ~ ., data = bfw_sig_snps[-1])
  bfw_partial <- lm(as.formula(paste("Bulliform_FW~", paste(snp_list[-i], collapse = "+"))), data = bfw_sig_snps[-1])
  temp_partial_R2 <- tibble(SNP = snp_list[i], 
                            Partial_R2 = round(partial.R2(bfw_partial, bfw_full), digits = 3))
  bfw_partial_R2 <- bind_rows(bfw_partial_R2, temp_partial_R2)
}
bfw_partial_R2 <- mutate(bfw_partial_R2, Trait = "BFW")

partial_R2s <- bind_rows(bff_partial_R2, bfw_partial_R2)

### % of variance explained by each SNP when all SNPs in model ###
### BFF ###
bff_aov_partial_r2 <- aov(Bulliform_FF ~ ., data = dplyr::select(bff_sig_snps, -Taxa)) %>%
  tidy() %>%
  mutate(partial_r2_all = sumsq/sum(sumsq), 
         Trait = "BFF") %>%
  dplyr::select(Trait, term, partial_r2_all) %>%
  dplyr::rename("SNP" = term)

### BFW ###
bfw_aov_partial_r2 <- aov(Bulliform_FW ~ ., data = dplyr::select(bfw_sig_snps, -Taxa)) %>%
  tidy() %>%
  mutate(partial_r2_all = sumsq/sum(sumsq), 
         Trait = "BFW") %>%
  dplyr::select(Trait, term, partial_r2_all) %>%
  dplyr::rename("SNP" = term)

aov_partial_r2 <- bind_rows(bff_aov_partial_r2, bfw_aov_partial_r2) %>%
  mutate(partial_r2_all = round(partial_r2_all, digits = 3))

### Create final table ###
sig_snps <- left_join(sig_snps, partial_R2s, by = c("Trait", "SNP")) %>%
  left_join(., aov_partial_r2, by = c("Trait", "SNP"))
