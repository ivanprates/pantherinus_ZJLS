###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To run genetic clustering analyses using sNMF;
### To make bar plots based on the ancestry coefficients.
### To make maps based on cluster assignments.

### PART 1: Getting ready ----

## Packages:
library(ggrepel)
library(LEA)
library(patchwork)
library(phytools)
library(scatterpie)
library(raster)
library(rnaturalearth)
library(sf)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_pantherinus/" ## Rhinellax.

## Working directory:
setwd(paste0(path, "sNMF"))
detach("package:here", unload = TRUE)
library(here)

## Creating directories to save results:
dir.create(path = here("data"))
dir.create(path = here("cross_entropy"))
dir.create(path = here("plots"))
dir.create(path = here("qmatrices"))

## Color palette for plots:
#palette <- c("#b6e2dc", "#0f85a0", "#5c4c35", "#c5352a", "#ffae57") ## Light green, cyano, brown, red, peach.
palette <- c("black", "#b6e2dc", "#0f85a0", "#c5352a", "#ffae57") ## Black, light green, cyano, red, peach.

## Subset:
palette <- palette[c(2:5)]
    
## PART 2: Filter genetic data ----

## Import genetic data:
dataset <- "pantherinus_R1_c85_n92_outfiles"
gendata <- read.table(paste0(path, "ipyrad/", dataset, "/usnps.012"), sep = "\t", row.names = 1, header = FALSE)
dim(gendata)

## Missing data per SNP:
calc_miss_SNP <- function(SNP) {
  
  ## Let's now estimate levels of missing data per sample:
  print(paste0("Now processing SNP ", SNP))
  
  ## How many sites have no data?
  miss <- table(gendata[, SNP] == 9)["TRUE"]
  
  ## What proportion this is from the total?
  if (is.na(miss) == TRUE) { 
    SNP_miss <- 0 
  } else { 
    SNP_miss <- round(miss/(dim(gendata)[1]), digits = 2) }
  
  ## Storing:
  SNP_miss <- data.frame(SNP = SNP, miss_prop = SNP_miss)
  SNP_miss$SNP <- as.character(SNP_miss$SNP)
  
  ## Return:
  return(SNP_miss)
  
} ## End of function.

## Apply:
SNP_df <- map_df(names(gendata), calc_miss_SNP)

## Save missing data information:
write.csv(SNP_df, file = "data/miss_data_per_SNP.csv", row.names = FALSE, quote = FALSE)

## Exclude SNPs based on maximum missing data:
miss_SNPs <- 0.3
SNPs_to_keep <- SNP_df[SNP_df$miss_prop <= miss_SNPs, ]

## Keep SNPs with less than max_miss missing data:
gendata <- gendata[names(gendata) %in% SNPs_to_keep$SNP]
dim(gendata)

## Sample IDs:
individuals <- read.table(paste0(path, "ipyrad/", dataset, "/usnps.012.indv"), header = FALSE)
gendata$SAMPLE_ID <- as.character(individuals$V1)

## Now, let's add locality information (including lat-longs) to the qmatrix:
sample_info <- read.csv(file = paste0(path, "sample_information/sample_information_pantherinus_ddRAD_2021-08-06.csv"), header = TRUE)

## Those in gendata:
sample_info <- sample_info[sample_info$SAMPLE_ID %in% gendata$SAMPLE_ID, ]

## How many samples per site?
site_list <- group_by(sample_info, LAT, LON) %>% group_split()
n_per_site <- purrr::map(.x = site_list, .f = nrow)

## Sites that have less or equal the maximum number of samples:
max_n <- 5
norm_sampled <- site_list[n_per_site <= max_n]
norm_sampled <- dplyr::bind_rows(norm_sampled)

## Oversampled sites:
over_sampled <- site_list[n_per_site > max_n]
over_sampled_filtered <- purrr::map_dfr(over_sampled, sample_n, size = max_n)

## Combine "normal" and filtered oversampled sites:
sample_info <- rbind(norm_sampled, over_sampled_filtered)

## Order:
sample_info <- arrange(sample_info, by = SAMPLE_ID)

## Keep only selected samples:
samples_keep <- intersect(gendata$SAMPLE_ID, sample_info$SAMPLE_ID)
gendata <- gendata[gendata$SAMPLE_ID %in% samples_keep, ]
dim(gendata)

## Remove ID info:
sample_IDs <- gendata$SAMPLE_ID
gendata <- gendata[setdiff(names(gendata), "SAMPLE_ID")]
dim(gendata)

## Save in geno format:
for (a in c(10, 50, 100, 200)) { write.geno(gendata, paste0("data/a", a, ".geno")) }

## PART 3: Function: Run sNMF ----

## Function:  
run_sNMF <- function(a, minK, maxK, project, rep) {
  
  ## Testing:
  #a <- 100 ; project <- "previous" ; rep <- 20 ; minK <- 3; maxK <- 6
  
  ## Running sNMF:
  if (project == "new") { project.snmf <- snmf(input.file = paste0("data/a", a, ".geno"), 
                          entropy = TRUE, ploidy = 2, project = "new", K = minK:maxK, alpha = a, repetitions = rep,
                          seed = 1984, CPU = 4) }
  
  ## Loading previously saved sNMF project:
  if (project == "previous") { project.snmf <- load.snmfProject(paste0("data/a", a, ".snmfProject")) }
  
  ## Showing summary of project results:
  snmf_summary <- summary(project.snmf)
  snmf_summary
  
  ## Selecting criterion to determine the best K:
  #crossEntropy <- snmf_summary$crossEntropy[1,] # Using min cross entropy across runs.
  crossEntropy <- snmf_summary$crossEntropy[2,] # Using mean cross entropy across runs.
  names(crossEntropy) <- gsub(x = names(crossEntropy), pattern = "K = ", replacement = "")
  
  ## Formatting:
  crossEntropy_df <- data.frame(K = names(crossEntropy), ce = crossEntropy)
  crossEntropy_df <- arrange(crossEntropy_df, by = ce)
  
  ## Selecting best K based on minimum (or mean) cross entropy among runs:
  bestK <- as.numeric(as.character(crossEntropy_df$K[which.min(crossEntropy_df$ce)]))
  K <- bestK
 
  ## Plotting mean entropy scores:
  png(filename = paste0("cross_entropy/a", a, "_K", K, ".png"))
    plot(y = crossEntropy, x = names(crossEntropy), cex = 1.2, col = "blue", pch = 19, xlab = paste0("K values"), ylab = "Cross-entropy")
    lines(y = crossEntropy, x = names(crossEntropy), col = "blue")
  dev.off() # Saving plot.
  
  ## Getting the cross entropy of all runs for K:
  ce <- cross.entropy(project.snmf, K = K)
    
  ## A few edits:
  ce_df <- as.data.frame(ce)
  ce_df$run <- 1:nrow(ce_df)
  ce_df <- ce_df[order(ce_df[, 1]), ]
  names(ce_df)[1] <- "ce"
    
  ## Selecting the run with the lowest cross-entropy for K = best K:
  bestrun <- which.min(ce)
    
  ## Getting the Q matrix:
  qmatrix <- as.data.frame(Q(project.snmf, K = K, run = bestrun))
      
  ## Replace column names:
  names(qmatrix) <- gsub(names(qmatrix), pattern = "V", replacement = "cluster_\\1")
  dim(qmatrix)
  
  ## Adding individual ID names:
  qmatrix$ID <- sample_IDs
    
  ## "Melt" dataframe using to assign samples to clusters based on max qscores:
  qmatrix_melt <- gather(qmatrix, key = assignment, value = coeff, 1:all_of(K))
    
  ## Assign specimens to cluster based on the highest qscore (coeff) values:
  cluster_assigned <- qmatrix_melt %>% group_by(ID) %>% top_n(n = 1, wt = coeff)
    
  ## Combine qmatrix and cluster assignments:
  qmatrix_c <- merge(qmatrix, cluster_assigned, by = "ID")
    
  ## Now, let's add locality information (including lat-longs) to the qmatrix:
  sample_info_all <- read.csv(file = paste0(path, "sample_information/sample_information_pantherinus_ddRAD_2021-08-06.csv"), header = TRUE)
  
  ## Let's only keep the individuals included in sNMF analyses.
  sample_info <- (sample_info_all[sample_info_all$SAMPLE_ID %in% qmatrix_c$ID, ])
  names(sample_info)[names(sample_info) == "SAMPLE_ID"] <- "ID" ## Renaming ID column for consistency across files.
    
  ## Merging:
  qmatrix_c <- merge(qmatrix_c, sample_info, by = "ID")
  
  ## Using a phylogeny to order samples in structure plots:
  tree <- read.nexus(paste0(path, "RAxML/RAxML_bipartitions_ddRAD_c85_2021-06-14.nex"))
  
  ## Obtaining the correct (i.e., rotated, ladderized) order of tips:
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  
  ## Keeping only columns of interest:
  assignments <- data.frame(ID = qmatrix_c$ID, assignment = qmatrix_c$assignment, coeff = qmatrix_c$coeff)
  
  ## Using tip labels to guide sample order in plots:
  phylo_order <- tree$tip.label[ordered_tips]
  assignments$ID <- factor(assignments$ID, levels = phylo_order)
  assignments <- arrange(assignments, ID)
  
  ## Saving new sample order: 
  assignments$plot_order <- 1:nrow(assignments)
  
  ## Keeping columns of interest:
  sample_order <- assignments[, c("ID", "plot_order")]
  
  ## Merge with qmatrix:
  qmatrix_c <- merge(qmatrix_c, sample_order, by = "ID")
  
  ## Save qmatrix:
  write.csv(qmatrix_c, file = here(paste0("qmatrices/a", a, ".csv")), row.names = FALSE)
  
} ## End of function.

## PART 4: Function: Plotting ancestry coefficients ----

## Creating function to plot with ggplot:
plot_sNMF_bars <- function(a) {

  ## Testing:
  #a <- 100
  
  ## Read corresponding qmatrix:
  qmatrix <- read.csv(file = paste0("qmatrices/a", a, ".csv"), header = TRUE)
  
  ## Simplify IDs:
  qmatrix$ID <- gsub(qmatrix$ID, pattern = "_Ct_pant", replacement = "")
  qmatrix$ID <- gsub(qmatrix$ID, pattern = "NA_([A-Z]+)([0-9]+)", replacement = "\\1 \\2")
  qmatrix$ID <- gsub(qmatrix$ID, pattern = "_", replacement = " ")
      
  ## Arrange:
  qmatrix <- arrange(qmatrix, desc(plot_order))

  ## Exclude admixed samples to set factor levels:
  qmatrix_l <- qmatrix[qmatrix$coeff > 0.95, ]

  ## Set factors levels to preserve order in plots:
  qmatrix$ID <- factor(qmatrix$ID, levels = qmatrix$ID)
    
  ## Number of clusters:
  K <- length(unique(qmatrix$assignment))
  
  ## Round qscores:
  for (column in paste0("cluster_", 1:K)) { qmatrix[column] <- round(qmatrix[column], 2) }
  
  ## "Melt":
  qmatrix_m <- gather(qmatrix, key = sNMF_cluster, value = qscores, 2:(K+1))
        
  ## Set factors levels to preserve order in plots:
  qmatrix_m$sNMF_cluster <- factor(qmatrix_m$sNMF_cluster, levels = unique(qmatrix_l$assignment))
    
  ## Creating stacked bar plot of ancestry coefficients:
  plot <- ggplot(data = qmatrix_m, aes(y = ID)) +
      
    ## Adding bars that represent ancestry coefficients:
    geom_bar(aes(x = qscores, fill = sNMF_cluster), 
             color = "white", ## Color of bar lines.
             size = 0.000001, ## Thickness of bar lines.
             stat = "identity", position = "fill", alpha = 0.85,
	           show.legend = FALSE) +
        
    ## Filling bars by cluster assigment:
    scale_fill_manual(values = rev(palette)) + 
          
    ## Adjusting labels:
    labs(x = "Ancestry proportions", y = "") +
          
    ## Adjusting limits of the x axis:
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = round(c(0, 1), 0)) +
            
    ## Changing theme:
    theme_minimal() +
  
    ## Theme parameters:
    theme(
          axis.text.y = element_text(color = "gray30", angle = 0, vjust = 0.5, hjust = 1, size = 8),
          #axis.text.y = element_blank(), ## Removing IDs from plot.
          axis.title.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid = element_blank(), 
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 0, margin = margin(0, 0, 0, 0)),
          plot.margin = margin(r = 15, l = 10, t = 0, b = 0),
          plot.title = element_blank())
	    
  ## Return:
  return(plot)
  
} ## End of function.

## PART 6: Function: Plotting maps ----

## Now creating a function to plot maps with ggplot:
plot_maps <- function(a) {
  
  ## Testing:
  #a <- 100
  
  ## Australia map:
  AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
  AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.
  
  ## Read qmatrix:
  qmatrix <- read.csv(file = here(paste0("qmatrices/a", a, ".csv")), header = TRUE)
       
  ## Arrange:
  qmatrix <- arrange(qmatrix, plot_order)

  ## Exclude admixed samples to set factor levels:
  qmatrix_l <- qmatrix[qmatrix$coeff >= 0.95, ]

  ## Set order of factors to preserve sample order in plots:
  qmatrix$assignment <- factor(qmatrix$assignment, levels = unique(qmatrix_l$assignment))
      
  ## Number of clusters:
  K <- length(unique(qmatrix$assignment))
  
  ## Change order of columns to keep cluster color order in map plots:
  qmatrix[c(paste0("cluster_", 1:K))] <- qmatrix[c(paste0("cluster_", 1:K))][levels(qmatrix$assignment)]
  
  ## Make qscores numeric:
  qmatrix[c(paste0("cluster_", 1:K))] <- apply(qmatrix[c(paste0("cluster_", 1:K))], 2, as.numeric)
  
  ## Calculating mean qscores per locality:
  qmatrix_mean <- qmatrix %>% group_by(assignment, LAT, LON) %>% summarise_at(c(paste0("cluster_", 1:K)), mean) 
  
  ## Round qscores:
  for (column in paste0("cluster_", 1:K)) { qmatrix_mean[column] <- round(qmatrix_mean[column], 2) }
  
  ## Plot map with ggplot:
  plot <- ggplot() +
        
    ## Adding baseline map:
    geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 0.5) +
        
    ## Pie chart of the ancestry coefficients: 
    geom_scatterpie(data = qmatrix_mean, aes(x = LON, y = LAT), 
                    pie_scale = 1, alpha = 0.85, color = "transparent",
                    cols = paste0("cluster_", 1:K)) +
                    
    ## Setting a bunch of aesthetic parameters:
    theme_void() +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_fill_manual(values = palette) +
        
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray40", size = 1),
          plot.title = element_blank(),
          plot.subtitle = element_blank())
          
  ## Return:
  #plot
  return(plot)
    
} ## End of function.

## PART 7: Run it all ----

## Function:
run_all <- function(a) {
  
  ## Run sNMF:
  run_sNMF(project = "new", a = a, minK = 3, maxK = 6, rep = 20)
  
  ## Create plots post-analysis:
  plot_b <- plot_sNMF_bars(a = a)
  plot_m <- plot_maps(a = a)

  ## Combining plots:
  plot_bm <- ( plot_b | plot_m ) + plot_layout(widths = c(1, 4))

  ## Saving plot:
  ggsave(plot = plot_bm, width = 34.5, height = 20, units = "cm", limitsize = FALSE, dpi = 500, filename = here(paste0("plots/a", a, "_", dim(gendata)[1], "_", dim(gendata)[2], ".jpg")))
  ggsave(plot = plot_bm, width = 34.5, height = 20, units = "cm", limitsize = FALSE, dpi = 500, filename = here(paste0("plots/a", a, "_", dim(gendata)[1], "_", dim(gendata)[2], ".pdf")))
  
} ## End of function.

## Run:
purrr::map(c(10, 50, 100, 200), run_all)

## End of script.
