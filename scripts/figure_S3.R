### Load libraries ###
##################################################
library(data.table)
library(tidyverse)
library(readr)

### Load data ###
##################################################
centromere_positions <- fread("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/Centromere_Positions.csv") %>%
  mutate(Chromosome = paste("Chr ", Chromosome, sep = "")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", 
                                                    "Chr 6", "Chr 7", "Chr 8", "Chr 9", "Chr 10")))

### Create functions ###
##################################################
ld_analysis_func <- function(path,
                             ld_threshold = 0.1, 
                             window_size = 1e6, 
                             step_size = 1e5){
  
  min_ld_decay <- fread(path) %>%
    filter(`R^2` <= ld_threshold) %>%
    group_by(Position1) %>%
    slice(which.min(Dist_bp))
  
  num_windows <- ceiling(((max(min_ld_decay$Position1) - min(min_ld_decay$Position1) + 1) - window_size)/step_size + 1)
  positions <- as.numeric()
  ld_avgs <- as.numeric()
  
  for(i in 1:num_windows){
    start <- min(min_ld_decay$Position1) + (i - 1)*step_size
    end <- start + window_size
    
    window_center <- (start + end)/2
    temp_avg <- mean(as.data.table(min_ld_decay)[Position1 >= start & Position1 <= end]$Dist_bp, na.rm = TRUE)
    
    positions <- c(positions, window_center)
    ld_avgs <- c(ld_avgs, temp_avg)
    
    print(paste(i, " of ", num_windows, " windows ", sep = ""))
  }
  return(tibble(Position = positions,
                Distance_bp = ld_avgs))
}

### Compute min physical distance for LD (R^2) to decay to 0.1 for each chromosome ###
##################################################
setwd("~/Box Sync/McNinch Projects/Bulliform_Cell_GWAS_Manuscript/data/ld_results")
path <- getwd()
file.names <- dir(path, pattern ="_final.txt")
sliding_ld_avgs <- tibble()
for(i in 1:length(file.names)){
  temp_sliding_ld_avgs <- ld_analysis_func(path = file.names[i]) %>%
    mutate(Chromosome = parse_number(file.names[i]))
  
  sliding_ld_avgs <- bind_rows(sliding_ld_avgs, temp_sliding_ld_avgs)
}
sliding_ld_avgs <- mutate(sliding_ld_avgs, Chromosome = paste("Chr ", Chromosome, sep = "")) %>%
  mutate(sliding_ld_avgs, Chromosome = factor(Chromosome, levels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", 
                                                                     "Chr 6", "Chr 7", "Chr 8", "Chr 9", "Chr 10")))

### Graph results ###
##################################################
ggplot(test, aes(x = Position/1000000, y = Distance_bp/1000)) +
  geom_point() +
  facet_wrap(vars(Chromosome), nrow = 5, ncol = 2, scales = "free") +
  geom_smooth(span = 0.05, se = FALSE, method = "loess") +
  geom_rect(data = centromere_positions, 
            aes(x = NULL, y = NULL, xmin = Start/1000000 - 1, xmax = End/1000000 + 1, ymin = -Inf, ymax = Inf), 
            alpha = 0.5, fill = "firebrick") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(size = 4), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 24), 
        strip.background = element_rect(fill = "white", size = 2, color = "black"), 
        strip.text = element_text(size = 24, face = "bold")) +
  labs(x = "Physical position (MB)", y = "Averge distance (Kb) between SNPs not in LD")