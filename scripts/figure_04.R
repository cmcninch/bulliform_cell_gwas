### Load libraries ###
##################################################
library(tidyverse)
library(data.table)
library(poppr)
library(ape)
library(lme4)
library(multcompView)
library(genetics)
library(GenomicRanges)
library(LDheatmap)
library(mapLD)
library(grid)
library(grDevices)

### Load data and format ###
##################################################
snp_positions <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/MAF_MissRate_Filtered_WiDiv_SNPs.map") %>%
  dplyr::select(V2, V1, V4)
colnames(snp_positions) <- c("SNP", "Chromosome", "Position")

all_snps <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/MAF_MissRate_Filtered_WiDiv_SNPs.xmat")

taxa <- all_snps$taxa

all_snps <- all_snps[, -1]

snp_info <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/snp_general_info.txt") %>%
  dplyr::select(`Site Name`, `Physical Position`, `Major Allele`, `Minor Allele`)
colnames(snp_info) <- c("SNP", "Position", "Major_Allele", "Minor_Allele")

Zm00001d007372_snps <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/Zm00001d007372_SNP_Positions.csv")

### Create functions ###
##################################################
blup_multiple_enviro_func <- function(df, trait, new_trait_name){
  blups <- ranef(lmer(get(trait) ~  (1 | SNP_Genotype_Name) + (1 | Location) + (1 | SNP_Genotype_Name:Location) + (1 | Rep:Location), data = df))$SNP_Genotype_Name %>%
    rownames_to_column(var = "Taxa") 
  colnames(blups) <- c("Taxa", new_trait_name)
  return(blups)
}

### Compute BLUPs ###
##################################################
### Trait enntry-means ###
bulliform_ff_entry_means <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/trait_data/Bulliform_Cell_Field_Counts.csv") %>%
  mutate(Bulliform_FF = Bulliform_Cell_Field_Count/1040 * 245) %>% # height of image is 1040 pixels; 245 pixels = 1mm
  group_by(SNP_Genotype_Name, Location, Rep, Row) %>%
  summarise(Bulliform_FF = mean(Bulliform_FF, na.rm = TRUE)) %>%
  drop_na()

### BLUPs ###
bff_both_blups <- blup_multiple_enviro_func(df = bulliform_ff_entry_means, trait = "Bulliform_FF", new_trait_name = "BFF_Both_Locations")

### LD at serine carboypeptidase1 ###
##################################################
snps <- Zm00001d007372_snps$SNP

### Compute LD ###
snp_calls <- cbind(data.frame(Taxa = taxa), 
                   as.data.frame(all_snps)[, which(colnames(all_snps) %in% snps)]) %>%
  gather("SNP", "SNP_Call", 2:ncol(.)) %>%
  left_join(., snp_info, by = "SNP") %>%
  mutate(Genotype = ifelse(SNP_Call == 0, paste(Major_Allele, Major_Allele, sep = "/"), 
                           ifelse(SNP_Call == 2, paste(Minor_Allele, Minor_Allele, sep = "/"), 
                                  paste(Major_Allele, Minor_Allele, sep = "/")))) %>%
  mutate(Allele1 = gsub('\\/.*', '', Genotype), 
         Allele2 = gsub('.*\\/', '', Genotype)) %>%
  mutate(Allele1 = ifelse(SNP_Call != 1, Allele1, NA), 
         Allele2 = ifelse(SNP_Call != 1, Allele2, NA))

temp_LD <- mapLD(SNPdata = as.data.frame(dplyr::select(snp_calls, Allele1, Allele2, Taxa, SNP)),
                 locusID.col = 'SNP',
                 subjectID.col = 'Taxa',
                 allele.cols = c('Allele1', 'Allele2'),
                 WhichGene = NA,
                 outgraph = NA)

### Graphing LD heatmaps ###
snp_calls <- dplyr::select(snp_calls, Taxa, SNP, Genotype) %>%
  spread(SNP, Genotype) %>%
  column_to_rownames(var = "Taxa") %>%
  mutate_all(., .funs = genotype)

MyHeatmap <- LDheatmap(snp_calls, 
                       as.numeric(gsub('.*\\_', '', colnames(snp_calls))), 
                       LDmeasure = "r",
                       title = NULL,
                       flip = TRUE,
                       add.map = FALSE,
                       SNP.name = snps,
                       color = c("#FF0000FF", "#FF4900FF", "#FF9200FF", "#FFDB00FF", "#FFFFE6FF"))

### Haplotypes at serine carboypeptidase1 ###
##################################################
snp_df <- as.data.frame(all_snps)[, c(1, which(colnames(all_snps) %in% c("rs2_230052852", "rs2_230053236")))] %>%
  mutate(Haplotype = paste(rs2_230052852, rs2_230053236, sep = "_")) %>%
  mutate(Hap_Code = ifelse(Haplotype == "0_0", 1, 
                           ifelse(Haplotype == "2_0", 2, 
                                  ifelse(Haplotype == "0_2", 3, 4))), 
         Taxa = taxa) %>%
  filter(Hap_Code != 4)
trait_snp_df <- left_join(snp_df, bff_both_blups, by = "Taxa")

aov_model <- aov(BFF_Both_Locations ~ Haplotype, data = trait_snp_df)

hsd <- TukeyHSD(aov_model, ordered = FALSE, conf.level = 0.95)

hsd_letters <- as.data.frame(multcompLetters(hsd$Haplotype[,4])[[1]]) %>%
  rownames_to_column(var = "Haplotype") %>%
  mutate(Y_Coord = c(0.7, 0.8, 0.5))
colnames(hsd_letters) <- c("Haplotype", "Letter", "Y_Coord")

### Graphing Haplotype bff blup distributions ###
ggplot(trait_snp_df, aes(x = Haplotype, y = BFF_Both_Locations)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", binwidth = 0.02, stackdir = "center", aes(fill = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_text(size = "24"), 
        axis.text = element_text(size = "20"), 
        panel.border = element_rect(size = 4), 
        legend.position = "none") +
  labs(y = "BFF BLUP") +
  geom_text(data = hsd_letters, 
            mapping = aes(label = Letter, 
                          x = Haplotype, 
                          y = Y_Coord), 
            size = 12)