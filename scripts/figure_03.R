### Load libraries and set working directory###
##################################################
library(tidyverse)
library(data.table)
library(lme4)
library(qvalue)
library(DESeq2)
library(GenomicRanges)
library(mapLD)
library(pheatmap)
library(RColorBrewer)

### Load data and format ###
##################################################
gff <- read_delim("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/Zea_mays.AGPv4.36.gff3", 
                  comment = "#", 
                  delim = "\t", 
                  progress = FALSE, 
                  col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  dplyr::select(X1, X4, X5, X9) %>%
  mutate(X9 = gsub('ID=gene:', '', X9)) %>%
  mutate(X9 = gsub(";.*", "", X9)) %>%
  arrange(X1, X4)
colnames(gff) <- c("Chromosome", "Gene_Start", "Gene_End", "Gene")

snp_positions <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/MAF_MissRate_Filtered_WiDiv_SNPs.map") %>%
  dplyr::select(V2, V1, V4)
colnames(snp_positions) <- c("SNP", "Chromosome", "Position")

all_snps <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/MAF_MissRate_Filtered_WiDiv_SNPs.xmat")

taxa <- all_snps$taxa

all_snps <- all_snps[, -1]

snp_info <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/snps/snp_general_info.txt") %>%
  dplyr::select(`Site Name`, `Physical Position`, `Major Allele`, `Minor Allele`)
colnames(snp_info) <- c("SNP", "Position", "Major_Allele", "Minor_Allele")

gene_ids <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/B73_V4_Gene_Name_Table.txt")

walley_library_info <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/Walley_2016_Library_Info.csv") %>%
  filter(Tissue %in% c("Leaf_8th_Mature", "Leaf_Growth_Zone", "Leaf_Stomatal_Division_Zone", "Leaf_Symmeterical_Division_Zone")) %>%
  mutate(Tissue = ifelse(Tissue == "Leaf_8th_Mature", "Mature_Leaf_FPKM", 
                         ifelse(Tissue == "Leaf_Growth_Zone", "Zone_3_FPKM",
                                ifelse(Tissue == "Leaf_Stomatal_Division_Zone", "Zone_2_FPKM",
                                       ifelse(Tissue == "Leaf_Symmeterical_Division_Zone", "Zone_1_FPKM", NA)))))

walley_fpkm <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/Walley_2016_FPKM_and_TPM_Maize_Atlas_B73_V4.csv") %>%
  left_join(gene_ids, ., by = c("Gene" = "GeneID")) %>%
  mutate(Library_Name = gsub(".tab", "", Library_Name)) %>%
  filter(Library_Name %in% walley_library_info$SRA_Run) %>%
  left_join(., walley_library_info, by = c("Library_Name" = "SRA_Run")) %>%
  group_by(GeneID, Tissue) %>%
  summarise(FPKM = mean(FPKM, na.rm = TRUE)) %>%
  spread(Tissue, FPKM)

count_data <- fread("~/Box Sync/McNinch Projects/Slice&Dice_RNA_Seq_Manuscript/data/Maize_Atlas/Walley_2016/Walley_2016_Gene_Count_Matrix_Maize_Atlas_V4_B73.csv") %>%
  dplyr::select(gene_id, walley_library_info$SRA_Run)

geneids <- fread("~/Box Sync/McNinch Projects/Slice&Dice_RNA_Seq_Manuscript/data/Maize_Atlas/Stepflug_2016/B73_V4_Gene_Name_Table.txt") %>%
  dplyr::rename("gene_id" = Gene)

count_data <- left_join(count_data, geneids, by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(GeneID, everything())

### Define functions ###
#################################################################################################################
# Differential expression analysis
DiffExp_func <- function(count_data, leaf_portion = c("Zone_1_FPKM", "Zone_2_FPKM", "Zone_3_FPKM", "Mature_Leaf_FPKM"), 
                         design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.34, FDR = 0.01, FC = 2){
  
  sample_info <- filter(walley_library_info, UQ(as.name(contrast[1])) %in% contrast[2:3], 
                        Tissue %in% leaf_portion) %>%
    arrange(UQ(as.name(contrast[1]))) %>% column_to_rownames(var = "SRA_Run")
  count_data <- dplyr::select(count_data, GeneID, rownames(sample_info)) %>% arrange(GeneID)
  max_low_expressed_libraries <- (1 - prop_libs_with_reads) * (ncol(count_data) -1)
  disqualified_genes <- count_data[rowSums(count_data < read_depth_cutoff) > max_low_expressed_libraries, ]$GeneID
  count_data <- count_data %>% column_to_rownames(var = "GeneID")
  ddsDF <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = design) 
  ddsDF <- DESeq(ddsDF)
  `%notIN%` <- Negate(`%in%`)
  res <- as.data.frame(results(ddsDF, contrast = contrast, alpha = FDR)) %>% rownames_to_column(var = "GeneID") %>% filter(padj <= FDR, GeneID %notIN% disqualified_genes)
  res <- bind_rows(filter(res, log2FoldChange >= log2(FC)), filter(res, log2FoldChange <= -log2(FC))) %>%
    dplyr::select(GeneID:log2FoldChange, padj) %>%
    mutate(Expression_Change = ifelse(log2FoldChange > 0, "Expression Decreases", "Expression Increases"))
  return(res)
}

### Compile significant SNPs based on q-value threshold ###
##################################################
setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Bulliform_GWAS_Results/")

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
  mutate(Trait = ifelse(Trait == "Bulliform_FF", "BFF",
                               ifelse(Trait == "Bulliform_FW", "BFW", "NA"))) %>%
  mutate(`-log10p.value` = round(`-log10p.value`, digits = 2), 
         `q_value` = formatC(`q_value`, format = "e", digits = 2))

### Define windows for each significant SNP to consider candidate genes ###
##################################################
window <- 20 ### Number of SNPs to the left or right of significant SNP
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

### Identify candidate genes ###
##################################################
windows <- left_join(sig_snps, 
                     dplyr::select(snps_in_LD, Significant_SNP, Upstream_Distance, Downstream_Distance), 
                     by = c("SNP" = "Significant_SNP")) %>%
  mutate(Upstream_Distance = ifelse(is.na(Upstream_Distance), 50000, 
                                    ifelse(Upstream_Distance <= 50000, 50000, Upstream_Distance)), 
         Downstream_Distance = ifelse(is.na(Downstream_Distance), 50000, 
                                      ifelse(Downstream_Distance <= 50000, 50000, Downstream_Distance))) %>%
  mutate(start = Position - Upstream_Distance, 
         end = Position + Downstream_Distance) %>%
  dplyr::select(Trait, SNP, Chromosome, Position, `-log10p.value`, q_value, start, end)

snp_ranges <- dplyr::select(sig_snps, Chromosome, Position, Trait, SNP) %>%
  left_join(., windows, by = c("Trait", "SNP", "Chromosome", "Position")) %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = start, end = end),
                  trait = Trait, 
                  SNP = SNP))

gene_ranges <- with(gff, GRanges(seqnames = Chromosome,
                                 ranges = IRanges(start = Gene_Start, end = Gene_End),
                                 ids = Gene))

snp_gene_overlap <- findOverlaps(gene_ranges, snp_ranges, ignore.strand = TRUE)

gene_indices <- snp_gene_overlap@from
snp_indices <- snp_gene_overlap@to

candidate_genes <- tibble(GeneID = gene_ranges@elementMetadata@listData$ids[gene_indices],
                          SNP = snp_ranges@elementMetadata@listData$SNP[snp_indices],
                          Trait = snp_ranges@elementMetadata@listData$trait[snp_indices]) %>%
  inner_join(., gff, by = c("GeneID" = "Gene")) %>%
  dplyr::select(Trait, SNP, GeneID, Gene_Start, Gene_End)

candidate_genes <- left_join(candidate_genes, sig_snps, by = c("Trait", "SNP")) %>%
  mutate(Gene_Start_Distance = Gene_Start - Position, 
         Gene_End_Distance = Gene_End - Position) %>%
  mutate(Gene_Distance_Kb = ifelse(abs(Gene_Start_Distance) < abs(Gene_End_Distance), 
                                   round(Gene_Start_Distance/1000, digits = 2), 
                                   round(Gene_End_Distance/1000, digits = 2))) %>%
  mutate(SNP_Within_Gene = ifelse(Position > Gene_Start & Position < Gene_End, "Yes", "No")) %>%
  mutate(Gene_Distance_Kb = ifelse(SNP_Within_Gene == "Yes", 0, Gene_Distance_Kb)) %>%
  dplyr::select(Trait, SNP, q_value, GeneID, Gene_Distance_Kb) 

### Differential expression analysis between leaf zones 1, 2 or 3 vs mature leaf ###
##################################################
zone1_degs <- DiffExp_func(count_data = count_data, leaf_portion = c("Zone_1_FPKM", "Mature_Leaf_FPKM"),
                           design = ~ Tissue, contrast = c("Tissue", "Zone_1_FPKM", "Mature_Leaf_FPKM")) %>%
  dplyr::select(GeneID, log2FoldChange, padj, Expression_Change)

zone2_degs <- DiffExp_func(count_data = count_data, leaf_portion = c("Zone_2_FPKM", "Mature_Leaf_FPKM"),
                           design = ~ Tissue, contrast = c("Tissue", "Zone_2_FPKM", "Mature_Leaf_FPKM")) %>%
  dplyr::select(GeneID, log2FoldChange, padj, Expression_Change)

zone3_degs <- DiffExp_func(count_data = count_data, leaf_portion = c("Zone_3_FPKM", "Mature_Leaf_FPKM"),
                           design = ~ Tissue, contrast = c("Tissue", "Zone_3_FPKM", "Mature_Leaf_FPKM")) %>%
  dplyr::select(GeneID, log2FoldChange, padj, Expression_Change)

# Compile summary
deg_summary <- tibble(GeneID = geneids$GeneID) %>%
  mutate(Zone_1 = ifelse(GeneID %in% zone1_degs$GeneID, 1, 0),
         Zone_2 = ifelse(GeneID %in% zone2_degs$GeneID, 1, 0),
         Zone_3 = ifelse(GeneID %in% zone3_degs$GeneID, 1, 0))

### Create final table ###
##################################################
candidate_genes <- left_join(candidate_genes, fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/GeneID_to_Locus_Lookup.csv"), by = "GeneID") %>%
  left_join(., dplyr::select(fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/gene_expression_data/B73v4_Gene_Annotations_Ensembl_58.csv"), GeneID, Gene_Description), by = "GeneID") %>%
  left_join(., walley_fpkm, by = "GeneID") %>%
  filter(Zone_1_FPKM > 1 |
           Zone_2_FPKM > 1 |
           Zone_3_FPKM > 1) %>%
  dplyr::select(Trait, SNP, GeneID, Gene_Distance_Kb, Gene_Description, Zone_1_FPKM, Zone_2_FPKM, Zone_3_FPKM, Mature_Leaf_FPKM) %>%
  mutate(Chromosome = as.numeric(gsub("rs", "", gsub('\\_.*', '', SNP))), 
         Position = as.numeric(gsub('.*\\_', '', SNP))) %>%
  dplyr::select(Trait, SNP, Chromosome, Position, everything()) %>%
  distinct(Trait, GeneID, .keep_all = TRUE) %>%
  mutate_at(.vars = vars(Zone_1_FPKM:Mature_Leaf_FPKM), round, digits = 2) %>%
  arrange(Trait, Chromosome, Position, Gene_Distance_Kb)

candidate_genes <- as.data.frame(candidate_genes)

bff_candidate_genes <- filter(candidate_genes, Trait == "BFF") %>%
  group_by(SNP) %>%
  mutate(SNP_Gene_Count = seq(n()))
bff_candidate_genes$Total_Gene_Count <- 1:nrow(bff_candidate_genes)

bff_candidate_gene_breaks <- sort((group_by(bff_candidate_genes, SNP) %>% dplyr::slice(base::which.max(SNP_Gene_Count)))$Total_Gene_Count)

### Graph Stats ###
#################################################################################################################
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)

colfunc <- colorRampPalette(c("yellow", "firebrick"))

# walley fpkm bff
pheatmap(bff_candidate_genes[, 7:10],
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         gaps_row = bff_candidate_gene_breaks,
         labels_row = bff_candidate_genes$SNP,
         scale = "row",
         cellwidth = 30,
         cellheight = 20,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         col = colfunc(15),
         legend = FALSE,
         legend_labels = NULL,
         border_color = "black", 
         na_col = "yellow", 
         display_numbers = FALSE)