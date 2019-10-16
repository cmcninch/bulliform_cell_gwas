### Load libraries and set working directory###
##################################################
library(tidyverse)
library(data.table)
library(lme4)
library(qvalue)
library(mapLD)

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

### Compile significant SNPs based on q-value threshold ###
##################################################
setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Flowering_Time_GWAS_Results/")

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

sig_snps <- dplyr::select(sig_snps, Trait, Chromosome, Position, SNP, `-log10p.value`, q_value) %>%
  mutate(`-log10p.value` = round(`-log10p.value`, digits = 2), 
         `q_value` = formatC(`q_value`, format = "e", digits = 2))

### Define windows for each significant SNP to determine LD-based haplotype blocks ###
##################################################
window <- 20 ### Number of SNPs to the left and right of significant SNP
snp_hits <- unique(sig_snps$SNP)

snps_in_LD <- tibble()
for(i in 1:length(snp_hits)){
  
  snp_of_interest_stats <- filter(snp_positions, SNP == snp_hits[i])
  snp_index <- which(snp_positions$SNP %in% snp_hits[i])
  neighboring_snps <- snp_positions[(snp_index - window):(snp_index + window),]
  neighboring_snps <- mutate(neighboring_snps, Significant_SNP = snp_hits[i])
  
  ### Compute LD blocks in windows ###
  snp_calls_in_window <- cbind(data.frame(Taxa = taxa), 
                               as.data.frame(all_snps)[, which(colnames(all_snps) %in% neighboring_snps$SNP)]) %>%
    gather("SNP", "SNP_Call", 2:ncol(.)) %>%
    left_join(., snp_info, by = "SNP") %>%
    mutate(Genotype = ifelse(SNP_Call == 0, paste(Major_Allele, Major_Allele, sep = "/"), 
                             ifelse(SNP_Call == 2, paste(Minor_Allele, Minor_Allele, sep = "/"), 
                                    paste(Major_Allele, Minor_Allele, sep = "/")))) %>%
    mutate(Allele1 = gsub('\\/.*', '', Genotype), 
           Allele2 = gsub('.*\\/', '', Genotype)) %>%
    mutate(Allele1 = ifelse(SNP_Call != 1, Allele1, NA), 
           Allele2 = ifelse(SNP_Call != 1, Allele2, NA))
  
  temp_LD <- mapLD(SNPdata = as.data.frame(dplyr::select(snp_calls_in_window, Allele1, Allele2, Taxa, SNP)),
                   locusID.col = 'SNP',
                   subjectID.col = 'Taxa',
                   allele.cols = c('Allele1', 'Allele2'),
                   WhichGene = NA,
                   outgraph = NA)
  
  snps_in_LD_index <- filter(temp_LD$LDblock, 
                             head <= filter(temp_LD$locusMap, LocusName == snp_hits[i])$LocusIndex[1] &
                               tail >= filter(temp_LD$locusMap, LocusName == snp_hits[i])$LocusIndex[1])
  
  if(nrow(snps_in_LD_index) == 0){
    temp_snps_in_LD <- tibble(Significant_SNP = snp_hits[i], 
                              SNPs_in_LD = NA)
  } else{
    temp_snps_in_LD <- tibble(Significant_SNP = snp_hits[i], 
                              SNPs_in_LD = as.character(temp_LD$locusMap[snps_in_LD_index$head[1]:snps_in_LD_index$tail[1], ]$LocusName))
  }
  
  temp_snps_in_LD <- mutate(temp_snps_in_LD, 
                            Sig_SNP_Position = as.numeric(gsub('.*\\_', '', Significant_SNP)), 
                            LD_SNP_Position = as.numeric(gsub('.*\\_', '', SNPs_in_LD)))
  
  temp_snps_in_LD <- mutate(temp_snps_in_LD, 
                            Upstream_Bound_Position = as.numeric(min(temp_snps_in_LD$LD_SNP_Position)), 
                            Downstream_Bound_Position = as.numeric(max(temp_snps_in_LD$LD_SNP_Position))) %>%
    mutate(Downstream_Bound_SNP = paste(gsub('\\_.*', '', Significant_SNP), "_", Downstream_Bound_Position, sep = ""), 
           Upstream_Bound_SNP = paste(gsub('\\_.*', '', Significant_SNP), "_", Upstream_Bound_Position, sep = "")) %>%
    mutate(Downstream_Distance = Downstream_Bound_Position - Sig_SNP_Position, 
           Upstream_Distance = Upstream_Bound_Position - Sig_SNP_Position) %>%
    dplyr::select(Significant_SNP, Upstream_Bound_SNP, Downstream_Bound_SNP, Upstream_Distance, Downstream_Distance) %>%
    distinct()
  
  snps_in_LD <- bind_rows(snps_in_LD, temp_snps_in_LD)
  print(i)
}