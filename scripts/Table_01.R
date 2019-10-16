### Load libraries ###
##################################################
library(tidyverse)
library(data.table)
library(lubridate)
library(lme4)
library(lmerTest)
library(broom.mixed)

### Create functions ###
##################################################
heritability_func <- function(df, mean, trait){
  
  model <- lmer(get(trait) ~  (1 | Genotype) + (1 | Location) + (1 | Genotype:Location) + (1 | Rep:Location), data = df) %>%
    tidy() %>%
    filter(is.na(std.error))
  
  var_g <- filter(model, group == "Genotype")$estimate[1]
  var_error <- filter(model, group == "Residual")$estimate[1]
  var_ge <- filter(model, group == "Genotype:Location")$estimate[1]
  
  return(var_g/(var_g + var_ge/2 + var_error/3*2))
}

repeatability_func <- function(df, trait){
  model <- lmer(get(trait) ~  (1 | Genotype) + (1 | Rep) + (1 | Genotype:Rep), data = df) %>%
    tidy() %>%
    filter(is.na(std.error))
  
  var_g <- filter(model, group == "Genotype")$estimate[1]
  var_error <- filter(model, group == "Residual")$estimate[1]
  
  return(var_g/(var_g + var_error/2))
}

### Load data and compute entry means ###
##################################################
bulliform_fw_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  group_by(Image, Genotype, SNP_Genotype_Name, Location, Row, Rep, Block) %>%
  summarise(Bulliform_FW = mean(Bulliform_Field_Area, na.rm = TRUE)/1388/950*1000) %>% # width of images is 1388 pixels; 950 pixels = 1mm
  group_by(Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_FW, na.rm = TRUE)) %>%
  drop_na()

bulliform_ff_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  group_by(Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FF = mean(Bulliform_Cell_Field_Count, na.rm = TRUE)/1040 * 245) %>% # height of images is 1040 pixels; 245 pixels = 1mm 
  drop_na()

flowering_times <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv") %>%
  mutate(Anthesis_Date = as.numeric(mdy(Anthesis_Date) - mdy("5/8/17")), 
         Silking_Date = as.numeric(mdy(Silking_Date) - mdy("5/8/17"))) %>%
  gather("Trait", "Score", 6:7) %>%
  group_by(Genotype, Rep, Row, Trait) %>%
  summarise(Score = mean(Score, na.rm = TRUE)) %>%
  spread(Trait, Score)

### Means and standard deviations ###
##################################################
bulliform_fw_stats <- tibble(trait = "bulliform_fw", 
                             mean = mean(bulliform_fw_entry_means$Bulliform_FW, na.rm = TRUE), 
                             sd = sd(bulliform_fw_entry_means$Bulliform_FW, na.rm = TRUE))

bulliform_ff_stats <- tibble(trait = "bulliform_ff", 
                             mean = mean(bulliform_ff_entry_means$Bulliform_FF, na.rm = TRUE), 
                             sd = sd(bulliform_ff_entry_means$Bulliform_FF, na.rm = TRUE))

dta_stats <- tibble(trait = "dta", 
                    mean = mean(flowering_times$Anthesis_Date, na.rm = TRUE),
                    sd = sd(flowering_times$Anthesis_Date, na.rm = TRUE))

dts_stats <- tibble(trait = "dts", 
                    mean = mean(flowering_times$Silking_Date, na.rm = TRUE),
                    sd = sd(flowering_times$Silking_Date, na.rm = TRUE))

combined_stats <- bind_rows(bulliform_fw_stats, bulliform_ff_stats, dta_stats, dts_stats)

### ANOVAs ###
##################################################
bff_multi_model <- aov(formula = Bulliform_FF ~  Genotype + Location + Genotype:Location + Location:Rep, data = bulliform_ff_entry_means)
summary(bff_multi_model)

bfw_multi_model <- aov(formula = Bulliform_FW ~  Genotype + Location + Genotype:Location + Location:Rep, data = bulliform_fw_entry_means)
summary(bfw_multi_model)

bff_single_model <- aov(Bulliform_FF ~  Genotype + Rep + Genotype:Rep, data = filter(bulliform_ff_entry_means, Location == "Ames"))
summary(bff_single_model)

bfw_single_model <- aov(Bulliform_FW ~  Genotype + Rep + Genotype:Rep, data = filter(bulliform_fw_entry_means, Location == "Ames"))
summary(bfw_single_model)

### Heritability and repeatability estimates ###
##################################################
heritability_func(df = bulliform_fw_entry_means, trait = "Bulliform_FW")
repeatability_func(df = filter(bulliform_fw_entry_means, Location == "Ames"), trait = "Bulliform_FW")
heritability_func(df = bulliform_ff_entry_means, trait = "Bulliform_FF")
repeatability_func(df = filter(bulliform_ff_entry_means, Location == "Ames"), trait = "Bulliform_FF")
repeatability_func(df = flowering_times, trait = "Anthesis_Date")
repeatability_func(df = flowering_times, trait = "Silking_Date")
