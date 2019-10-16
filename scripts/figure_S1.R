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
library(GGally)
library(ggcorrplot)

### Create functions ###
##################################################
blup_multiple_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | Genotype) + (1 | Location) + (1 | Genotype:Location) + (1 | Rep:Location), data = df))$Genotype %>%
    rownames_to_column(var = "Genotype") 
  colnames(blups) <- c("Genotype", trait)
  return(blups)
}

blup_single_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | Genotype) + (1 | Rep) + (1 | Genotype:Rep), data = df))$Genotype %>%
    rownames_to_column(var = "Genotype")
    colnames(blups) <- c("Genotype", trait)
  return(blups)
}

### Load data and compute blups ###
##################################################
### Trait entry-means ###
bulliform_fw_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  group_by(Image, Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_Field_Area, na.rm = TRUE)/1388/950*1000) %>% # width of images is 1388 pixels; 950 pixels = 1mm
  group_by(Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_FW, na.rm = TRUE))

bulliform_ff_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  mutate(Bulliform_FF = Bulliform_Cell_Field_Count/1040 * 245) %>% # height of images is 1040 pixels; 245 pixels = 1mm
  ungroup() %>%
  group_by(Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FF = mean(Bulliform_FF, na.rm = TRUE))

flowering_times_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv") %>%
  mutate(DTA = as.numeric(mdy(Anthesis_Date) - mdy("5/8/17")), 
         DTS = as.numeric(mdy(Silking_Date) - mdy("5/8/17"))) %>%
  dplyr::select(Genotype, Row, Rep, Block, DTA, DTS) %>%
  gather("Trait", "Score", 5:6) %>%
  group_by(Genotype, Row, Rep, Block, Trait) %>%
  summarise(Score = mean(Score, na.rm = TRUE)) %>%
  spread(Trait, Score)

### BLUPs ###
bulliform_ff_blups <- blup_multiple_enviro_func(df = bulliform_ff_entry_means, trait = "Bulliform_FF")
bulliform_fw_blups <- blup_multiple_enviro_func(df = bulliform_fw_entry_means, trait = "Bulliform_FW")
dta_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTA")
dts_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTS")

combined_blups <- left_join(bulliform_ff_blups, bulliform_fw_blups, by = "Genotype") %>%
  left_join(., dta_blups, by = "Genotype") %>%
  left_join(., dts_blups, by = "Genotype")

colnames(combined_blups) <- c("Genotype", "BFF", "BFW", "DTA", "DTS")

### Graph stats ###
##################################################
### Pairwise correlations ###
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black") +
    geom_smooth(method = method, color = "red", ...)
  p
}

ggpairs(combined_blups[-1], lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "blue")),
        upper = list(continuous = wrap("cor", size = 12))) + 
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 24), 
        panel.border = element_rect(size = 2), 
        axis.text = element_text(size = 12))

### Significance of correlations ###
cor_pmat(combined_blups[, -1], use = "complete.obs")
