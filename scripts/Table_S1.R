### Load libraries ###
##################################################
library(tidyverse)
library(data.table)
library(stats)
library(broom)
library(scales)
library(lubridate)
library(lme4)
library(broom)

### Create functions ###
##################################################
blup_multiple_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Location) + (1 | SNP_Genotype_Name:Location) + (1 | Rep:Location), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "SNP_Genotype_Name") 
  colnames(blups) <- c("SNP_Genotype_Name", trait)
  return(blups)
}

blup_single_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Rep) + (1 | SNP_Genotype_Name:Rep), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "SNP_Genotype_Name")
    colnames(blups) <- c("SNP_Genotype_Name", trait)
  return(blups)
}

### Load data and compute blups ###
##################################################
### Summary of where traits measured ###
bulliform_ff_summary <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  distinct(Genotype, Location) %>%
  mutate(Present = "Yes") %>%
  spread(Location, Present)
colnames(bulliform_ff_summary) <- c("Genotype", "BFF_in_Ames", "BFF_in_Madison")

bulliform_fw_summary <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  distinct(Genotype, Location) %>%
  mutate(Present = "Yes") %>%
  spread(Location, Present)
colnames(bulliform_fw_summary) <- c("Genotype", "BFW_in_Ames", "BFW_in_Madison")

genotype_lookup <- distinct(bind_rows(distinct(fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv"), Genotype, SNP_Genotype_Name),
                             distinct(fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv"), Genotype, SNP_Genotype_Name),
                             distinct(fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv"), Genotype, SNP_Genotype_Name)),
                            Genotype, SNP_Genotype_Name) %>%
  mutate(SNP_Genotype_Name = ifelse(is.na(SNP_Genotype_Name), "No SNPs", SNP_Genotype_Name))

### Trait entry-means ###
bulliform_fw_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  group_by(Image, SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_Field_Area, na.rm = TRUE)/1388/950*1000) %>% # width of images is 1388 pixels; 950 pixels = 1mm
  group_by(SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_FW, na.rm = TRUE))

bulliform_ff_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  mutate(Bulliform_FF = Bulliform_Cell_Field_Count/1040 * 245) %>% # height of images is 1040 pixels; 245 pixels = 1mm
  ungroup() %>%
  group_by(SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FF = mean(Bulliform_FF, na.rm = TRUE))

flowering_times_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv") %>%
  mutate(DTA = as.numeric(mdy(Anthesis_Date) - mdy("5/8/17")), 
         DTS = as.numeric(mdy(Silking_Date) - mdy("5/8/17"))) %>%
  dplyr::select(SNP_Genotype_Name, Row, Rep, Block, DTA, DTS) %>%
  gather("Trait", "Score", 5:6) %>%
  group_by(SNP_Genotype_Name, Row, Rep, Block, Trait) %>%
  summarise(Score = mean(Score, na.rm = TRUE)) %>%
  spread(Trait, Score)

### BLUPs ###
bulliform_ff_blups <- blup_multiple_enviro_func(df = bulliform_ff_entry_means, trait = "Bulliform_FF")
bulliform_fw_blups <- blup_multiple_enviro_func(df = bulliform_fw_entry_means, trait = "Bulliform_FW")
dta_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTA")
dts_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTS")

combined_blups <- left_join(bulliform_ff_blups, bulliform_fw_blups, by = "SNP_Genotype_Name") %>%
  left_join(., dta_blups, by = "SNP_Genotype_Name") %>%
  left_join(., dts_blups, by = "SNP_Genotype_Name")

colnames(combined_blups) <- c("SNP_Genotype_Name", "BFF_BLUP", "BFW_BLUP", "DTA_BLUP", "DTS_BLUP")

### Combine sumarries ###
##################################################
trait_measurement_summary <- full_join(genotype_lookup, bulliform_ff_summary, by = "Genotype") %>%
  full_join(., bulliform_fw_summary, by = "Genotype") %>%
  full_join(., combined_blups, by = "SNP_Genotype_Name") %>%
  mutate(DTA_in_Ames = ifelse(!is.na(DTA_BLUP),  "Yes", "No"),
         DTS_in_Ames = ifelse(!is.na(DTS_BLUP),  "Yes", "No"))
