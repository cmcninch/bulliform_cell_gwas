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
density_plot_func <- function(trait, legend_positions, x_axis_title, x_scales){
  df <- filter(combined_traits, Trait == trait) %>%
    dplyr::select(-Trait) %>%
    gather("Location", "Score", 2:3)

  ggplot(df, aes(x = Score, fill = Location)) +
    geom_density(size = 2, alpha = 0.4) +
    theme_bw() +
    scale_x_continuous(limits = c(min(x_scales), max(x_scales)), breaks = x_scales) +
    theme(axis.text = element_text(size = 24), 
          axis.title = element_text(siz = 30),
          panel.border = element_rect(size = 4), 
          legend.position = legend_positions, 
          legend.background = element_rect(color = "black", size = 2), 
          legend.text = element_text(size = 24)) +
    labs(x = x_axis_title, 
         fill = NULL)
}

### Load data ###
##################################################
bulliform_ff <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  mutate(Bulliform_FF = Bulliform_Cell_Field_Count/1040 * 245) %>%# height of images is 1040 pixels; 245 pixels = 1mm
  ungroup() %>%
  dplyr::select(-SNP_Genotype_Name)

bulliform_fw <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Areas.csv") %>%
  group_by(Image, Genotype, Location, Rep, Row) %>%
  summarise(Bulliform_FW = mean(Bulliform_Field_Area, na.rm = TRUE)/1388/950*1000) %>%# width of images is 1388 pixels; 950 pixels = 1mm
  ungroup()

### Compute stats ###
##################################################
### Genotype-Location averages ###
bulliform_ff_avgs <- group_by(bulliform_ff, Genotype, Location) %>%
  summarise(Bulliform_FF = mean(Bulliform_FF, na.rm = TRUE))

bulliform_fw_avgs <- group_by(bulliform_fw, Genotype, Location) %>%
  summarise(Bulliform_FW = mean(Bulliform_FW, na.rm = TRUE))

combined_traits <- bind_rows(spread(bulliform_ff_avgs, Location, Bulliform_FF) %>% mutate(Trait = "BFF"),
                             spread(bulliform_fw_avgs, Location, Bulliform_FW) %>% mutate(Trait = "BFW"))

### Graph stats ###
##################################################
### Density plots ###
density_plot_func(trait = "BFF", legend_positions = c(0.2, 0.8), x_axis_title = "BFF (fields per mm)", x_scales = c(1, 2, 3))

density_plot_func(trait = "BFW", legend_positions = c(0.175, 0.8), x_axis_title = "BFW (um)", x_scales = c(30, 45, 60, 75))
