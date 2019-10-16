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
blup_multiple_enviro_func <- function(df, trait){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Location) + (1 | SNP_Genotype_Name:Location) + (1 | Rep:Location), data = df))$SNP_Genotype_Name %>%
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

combined_blups <- full_join(bulliform_ff_blups, bulliform_fw_blups, by = "Taxa")
combined_blups[combined_blups == "NaN"] <- NA

combined_blups <- left_join(tibble(Taxa = desc@description$rowNames), combined_blups, by = "Taxa")

### Check all files have taxa and snps ordered exactly the same ###
identical(myCV$taxa, combined_blups$Taxa)
identical(myGM$SNP, desc@description$colNames)

### Run GWAS ###
##################################################
dir.create("../Bulliform_GWAS_Results")
setwd("../Bulliform_GWAS_Results")

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

setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Bulliform_GWAS_Results/")

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
  mutate(Trait = ifelse(Trait == "Bulliform_FF", "BFF",
                        ifelse(Trait == "Bulliform_FW", "BFW", "NA"))) %>%
  mutate(Adjusted_Position = Position + Position_Adjust, 
         Color = ifelse(`q-value` <= 0.01, "firebrick1", 
                        ifelse(`q-value` > 0.01 & `q-value` <= 0.05, "darkorange",
                               ifelse(Chromosome %in% c(1, 3, 5, 7, 9), "blue4", "blue1")))) %>%
  dplyr::select(Trait, SNP, Adjusted_Position, p.value, `-log10(p-value)`, `q-value`, Color, Breaks)

### Make manhattan plots ###
tiff("~/BFF_manhattan.tiff", res = 300, width = 6, height = 4, units = "in")
myplot <- manhattan_plot_func(df = snp_values, trait = "BFF", x_label = NULL,
                              y_lims = c(0, 11), y_breaks = c(0, 2, 4, 6, 8, 10), y_labs = c("0", "2", "4", "6", "8", "10"))
print(myplot)
dev.off()

tiff("~/BFW_manhattan.tiff", res = 300, width = 6, height = 4, units = "in")
myplot <- manhattan_plot_func(df = snp_values, trait = "BFW", x_label = NULL,
                              y_lims = c(0, 11), y_breaks = c(0, 2, 4, 6, 8, 10), y_labs = c("0", "2", "4", "6", "8", "10"))
print(myplot)
dev.off()

### Make qq plots ###
tiff("~/BFF_qq.tiff", res = 300, width = 4, height = 4, units = "in")
par(mar = c(5,5,1,1))
myplot <- qqPlot(filter(snp_values, Trait == "BFF")$p.value, truncate = FALSE, main = NULL, cex.lab = 1.5, cex.axis = 1.25, cex.main = 2)
print(myplot)
dev.off()

tiff("~/BFW_qq.tiff", res = 300, width = 4, height = 4, units = "in")
par(mar = c(5,5,1,1))
myplot <- qqPlot(filter(snp_values, Trait == "BFW")$p.value, truncate = FALSE, main = NULL, cex.lab = 1.5, cex.axis = 1.25, cex.main = 2)
print(myplot)
dev.off()