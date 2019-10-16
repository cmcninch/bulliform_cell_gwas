### Load libraries and set working directory###
##################################################
library(tidyverse)
library(lubridate)
library(data.table)
library(lme4)
library(broom)
library(FarmCPUpp)
library(bigmemory)
library(qvalue)
library(GWASTools)
library(grDevices)

setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps")

### Create functions ###
##################################################
blup_single_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Rep) + (1 | SNP_Genotype_Name:Rep), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "Taxa")
  colnames(blups) <- c("Taxa", trait)
  return(blups)
}

### Make manhattan plots ###
manhattan_plot_func <- function(df, trait, y_lims, y_breaks, y_labs, x_label){
  ggplot(filter(df, Trait == trait), 
         aes(x = Adjusted_Position, y = `-log10(p-value)`, color = Color, size = Color)) +
    geom_point() +
    guides(color = FALSE, size = FALSE) +
    scale_x_continuous(breaks = chromosome_sizes$Breaks, labels = 1:10, 
                       limits = c(-1000, 2106339117), expand = expand_scale(mult = c(0))) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_labs) +
    scale_color_manual(values = c("blue4", "blue1", "darkorange", "firebrick1")) +
    scale_size_manual(values = c(1, 1, 3, 3)) +
    theme_bw() +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size = 24),
          panel.grid = element_blank(), 
          panel.border = element_rect(size = 2)) +
    labs(x = x_label, 
         y = "-log10(p-value)")
}

### Load data ###
##################################################
### SNP data ###
myGM <- as.data.frame(fread("./MAF_MissRate_Filtered_WiDiv_SNPs.map")) %>%
  dplyr::select(V2) %>%
  mutate(Chromosome = as.numeric(gsub("rs", "", gsub("\\_.*", "", V2))), 
         Position = as.numeric(gsub('.*\\_', '', V2)))
colnames(myGM) <- c("SNP", "Chromosome", "Position")

myGD <- read.big.matrix("./MAF_MissRate_Filtered_WiDiv_SNPs.xmat",
                        sep = "\t",
                        header = TRUE,
                        col.names = myGM$SNP, 
                        has.row.names = TRUE,
                        ignore.row.names = FALSE,
                        type = "double",
                        backingfile = "MAF_MissRate_Filtered_WiDiv_SNPs",
                        descriptorfile = "MAF_MissRate_Filtered_WiDiv_SNPs.desc")

dput(describe(myGD), "MAF_MissRate_Filtered_WiDiv_SNPs.desc")

desc <- dget("MAF_MissRate_Filtered_WiDiv_SNPs.desc")
myGD <- attach.big.matrix(desc)

### Covariate matrices ###
myCV <- read_csv("GAPIT.PCA.csv") # Q matrix derived with default parameters of the GAPIT function

### Phenotype data ###
flowering_times_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Flowering_Traits.csv") %>%
  mutate(DTA = as.numeric(mdy(Anthesis_Date) - mdy("5/8/17")), 
         DTS = as.numeric(mdy(Silking_Date) - mdy("5/8/17"))) %>%
  dplyr::select(Genotype, Row, Rep, Block, DTA, DTS) %>%
  gather("Trait", "Score", 5:6) %>%
  group_by(Genotype, Row, Rep, Block, Trait) %>%
  summarise(Score = mean(Score, na.rm = TRUE)) %>%
  spread(Trait, Score)

### BLUPs ###
dta_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTA")
dts_blups <- blup_single_enviro_func(df = flowering_times_entry_means, trait = "DTS")

combined_blups <- full_join(dta_blups, dts_blups, by = "Taxa")
combined_blups <- left_join(tibble(Taxa = myCV$taxa), combined_blups, by = "Taxa")

### Check all files have taxa and snps ordered exactly the same ###
identical(myCV$taxa, combined_blups$Taxa)
identical(myGM$SNP, desc@description$colNames)

### Run GWAS ###
##################################################
dir.create("../Flowering_Time_GWAS_Results")
setwd("../Flowering_Time_GWAS_Results")

for(i in 2:ncol(combined_blups)){
  temp <- farmcpu(Y = combined_blups[, c(1, i)],
                  GD = myGD,
                  GM = myGM,
                  CV = dplyr::select(myCV, PC1:PC3),
                  method.bin = "optimum")
write_results(temp)
}

### Prepare SNPs ###
##################################################
chromosome_sizes <- tibble(Chromosome = 1:10, 
                           Size = c(307041717, 244442276, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 150982314)) %>%
  mutate(Position_Adjust = cumsum(Size)) %>%
  mutate(Position_Adjust = Position_Adjust - Size) %>%
  mutate(Breaks = Size/2 + Position_Adjust)

setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Flowering_Time_GWAS_Results")

path <- getwd()
out.file <- tibble()
file.names <- dir(path, pattern =".GWAS.csv")
for(i in 1:length(file.names)){
  file <- read_csv(file.names[i], 
                   col_names = TRUE) %>%
    mutate(Trait = gsub("\\.GWAS\\.csv", "", file.names[i]), 
           `-log10(p-value)` = as.numeric(-log10(p.value)))
  file$`q-value` <- qvalue(file$p.value)$qvalues
  file <- dplyr::select(file, Trait, SNP, Chromosome, Position, p.value, `-log10(p-value)`, `q-value`)
  out.file <- bind_rows(out.file, file)
}
snp_values <- left_join(out.file, chromosome_sizes, by = "Chromosome") %>%
  mutate(Adjusted_Position = Position + Position_Adjust, 
         Color = ifelse(`q-value` <= 0.01, "firebrick1", 
                        ifelse(`q-value` > 0.01 & `q-value` <= 0.05, "darkorange",
                               ifelse(Chromosome %in% c(1, 3, 5, 7, 9), "blue4", "blue1")))) %>%
  dplyr::select(Trait, SNP, Adjusted_Position, p.value, `-log10(p-value)`, `q-value`, Color, Breaks)

### Make manhattan plots ###
tiff("~/DTA_manhattan.tiff", res = 300, width = 6, height = 4, units = "in")
myplot <- manhattan_plot_func(df = snp_values, trait = "DTA", x_label = NULL,
                              y_lims = c(0, 13), y_breaks = c(0,  4, 8, 12), y_labs = c("0", "4", "8", "12"))
print(myplot)
dev.off()

tiff("~/DTS_manhattan.tiff", res = 300, width = 6, height = 4, units = "in")
myplot <- manhattan_plot_func(df = snp_values, trait = "DTS", x_label = NULL,
                              y_lims = c(0, 17), y_breaks = c(0, 4, 8, 12, 16), y_labs = c("0", "4", "8", "12", "16"))
print(myplot)
dev.off()

### Make qq plots ###
tiff("~/DTA_qq.tiff", res = 300, width = 4, height = 4, units = "in")
par(mar = c(5,5,1,1))
myplot <- qqPlot(filter(snp_values, Trait == "DTA")$p.value, truncate = FALSE, main = NULL, cex.lab = 1.5, cex.axis = 1.25, cex.main = 2)
print(myplot)
dev.off()

tiff("~/DTS_qq.tiff", res = 300, width = 4, height = 4, units = "in")
par(mar = c(5,5,1,1))
myplot <- qqPlot(filter(snp_values, Trait == "DTS")$p.value, truncate = FALSE, main = NULL, cex.lab = 1.5, cex.axis = 1.25, cex.main = 2)
print(myplot)
dev.off()
