###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To perform genetic PCA based on the ddRAD data.
### To perform NMDS based on morphological characters scored in museum subspecies of Ctenotus pantherinus.
### To plot the results.

## PART 1: Getting ready ----

## Packages:
library(ggforce)
library(LEA)
library(patchwork)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_pantherinus/")
detach("package:here", unload = TRUE)
library(here)

## Shapes:
shapes <- c(22, 24, 21, 23)

## Color palette:
palette <- c("#ae5f23", "#f3e2c6", "#15b6ad", "#152868")
palette <- palette[c(1, 3, 2, 4)]

## PART 2: Genetic PCA ----

## Process SNP data into geno format:
gendata <- read.table("VCFtools/snps.012", sep = "\t", row.names = 1)
#write.geno(gendata, output.file = "VCFtools/snps.geno")

## Location of geno file:
geno <- here("VCFtools/snps.geno")

## Run PCA analysis:
PCA <- pca(geno)

## Displaying information on the PCA analysis:
summary(PCA)

## Plotting eigenvalues:
plot(PCA, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")

## Storing the proportion of variation represented by selected PC axes"
PC_perc <- vector("list", 4) # We'll use the four first PC axes and save them into this list.
for (p in c(1:4)) { PC_perc[[p]] <- round((summary(PCA)[2, p]*100), digits = 0) 
names(PC_perc) <- c(paste0("PC", rep(1:4))) }

## Saving genetic PCA axes to plot later:
pcadata <- as.data.frame(PCA$projections)
pcadata <- pcadata[, 1:4]
names(pcadata) <- c(paste0("PC", rep(1:4)))

## Add sample IDs:
SAMPLE_ID <- data.frame(SAMPLE_ID = read.table("VCFtools/snps.012.indv", header = FALSE))
names(SAMPLE_ID) <- "SAMPLE_ID"
pcadata$SAMPLE_ID <- SAMPLE_ID$SAMPLE_ID

## ddRAD clade and subspecies assignments:
clade_df <- read.csv(file = "data/assignments_ddRAD_clades.csv", header = TRUE)
sample_info <- read.csv(file = "sample_information/sample_information_pantherinus_ddRAD_2021-08-06.csv", header = TRUE)
clade_ll <- merge(sample_info, clade_df, by = "SAMPLE_ID")

## Merge pca data with subspecies assignments:
pcadata <- merge(pcadata, clade_ll, by = "SAMPLE_ID")
names(pcadata)[names(pcadata) == "SUBSPECIES"] <- "Subspecies"

## Levels:
pcadata$clade <- factor(pcadata$clade, levels = sort(unique(pcadata$clade)))
pcadata$subspecies <- factor(pcadata$Subspecies, levels = sort(unique(pcadata$Subspecies)))

## Plotting:
plot_PCA <- ggplot(pcadata) +
  
  ## Plotting ellipses for each group:
  geom_mark_ellipse(data = filter(pcadata, clade == "ddRAD clade 1"), aes(x = PC1, y = PC2), color = "gray70", expand = unit(3, "mm"), size = 0.3, alpha = 0.05) +
  geom_mark_ellipse(data = filter(pcadata, clade == "ddRAD clade 2"), aes(x = PC1, y = PC2), color = "gray70", expand = unit(3, "mm"), size = 0.3, alpha = 0.05) +
  geom_mark_ellipse(data = filter(pcadata, clade == "ddRAD clade 3"), aes(x = PC1, y = PC2), color = "gray70", expand = unit(3, "mm"), size = 0.3, alpha = 0.05) +
  geom_mark_ellipse(data = filter(pcadata, clade == "ddRAD clade 4"), aes(x = PC1, y = PC2), color = "gray70", expand = unit(3, "mm"), size = 0.3, alpha = 0.05) +
  
  ## Labels:
  geom_text(aes(label = "ddRAD clade 1", y = -50, x = -50), color = "gray70", size = 4) +
  geom_text(aes(label = "ddRAD clade 2", y = 14, x = -48), color = "gray70", size = 4) +
  geom_text(aes(label = "ddRAD clade 3", y = 45, x = 11), color = "gray70", size = 4) +
  geom_text(aes(label = "ddRAD clade 4", y = -26, x = 14), color = "gray70", size = 4) +
  
  ## Configuring point size and shape:
  geom_point(aes(x = PC1, y = PC2, fill = .data[["Subspecies"]], shape = .data[["Subspecies"]]), 
             size = 3, color = "black", alpha = 0.75, show.legend = FALSE) +
  
  ## Labels:
  #geom_text(aes(x = PC1, y = PC2, label = .data[[group_var]])) +
  
  ## Expand plot limits:
  expand_limits(y = c(min(pcadata$PC2)-2, max(pcadata$PC2)+2), x = c(min(pcadata$PC1)-2, max(pcadata$PC1)+2)) +
  
  ## Axis labels:
  labs(y = paste0("PC2 (", PC_perc$PC2, "%)"), 
       x = paste0("PC1 (", PC_perc$PC1, "%)")) +
  
  ## Setting color scheme and, if needed, removing legend:
  scale_fill_manual(values = palette) + 
  scale_color_manual(values = palette) +
  scale_shape_manual(values = c(21:24)) +
  
  ## Also setting values on axes:
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  
  ## Other edits:
  theme_bw() +
  theme(plot.margin = margin(l = 5, r = 20),
        panel.border = element_rect(size = 1, colour = "gray20"),
        axis.text = element_text(color = "gray20", size = 12),
        axis.title = element_text(color = "gray20", size = 14),
        panel.grid = element_blank()) ## Removing background grid.
        
## Check:
plot_PCA

## Saving plot:
#ggsave(filename = paste0("figures/gen_PCA.jpg"), plot = plot_PCA, width = 16, height = 12, dpi = 500, units = "cm")
#ggsave(filename = paste0("figures/gen_PCA.pdf"), plot = plot_PCA, width = 16, height = 12, dpi = 500, units = "cm")

## Remove pca-related files:
unlink("*.pcaProject")
unlink("ipyrad/pantherinus_R1_c85_n92_outfiles/*.lfmm")
unlink("ipyrad/pantherinus_R1_c85_n92_outfiles/*.pca", recursive = TRUE)

## PART 3: NMDS ----

## Data:
char_df <- read.csv(file = "morphology/pantherinus_morphological_data_organized_ssps_added.csv", header = TRUE)

## Minor edit to prefrontal:
char_df$Prefrontal <- as.character(char_df$Prefrontal)
char_df$Prefrontal[char_df$Prefrontal == "narrowly separated" & !is.na(char_df$Prefrontal)] <- "narrow"
char_df$Prefrontal <- as.factor(char_df$Prefrontal)

morph_df <- char_df[c("REGNO", "Taxon", 
                      "No_midbody_rows", "No_lamellae", "Snout_vent_length", "SupLab",
                      "Dorsal_pattern", "Leg_pattern", "Digital_lamellae", "Plantar_scales", "Nasals", "Prefrontal")]
morph_df[morph_df == ""] <- NA

## Discrete variables:
disc_df <- morph_df[c("REGNO", "Dorsal_pattern", "Leg_pattern", "Digital_lamellae", "Plantar_scales")] ## "Nasals", "Prefrontal" not variable.

## 1. If ordering discrete character states (as per Rabosky et al. 2014 MPE):

## Order levels:
disc_wide <- disc_df
disc_wide$Dorsal_pattern <- factor(disc_wide$Dorsal_pattern, levels = c("Vert. hiatus", "Typical", "Wide border", "Vert. stripe", "Lines"), labels = c(1:length(setdiff(unique(disc_wide$Dorsal_pattern), NA))))
disc_wide$Digital_lamellae <- factor(disc_wide$Digital_lamellae, levels = c("Single broad", "Single fine", "Triple fine"), labels = c(1:length(setdiff(unique(disc_wide$Digital_lamellae), NA))))
disc_wide$Leg_pattern <- factor(disc_wide$Leg_pattern, levels = c("Spots", "Combo", "Lines"), labels = c(1:length(setdiff(unique(disc_wide$Leg_pattern), NA))))
disc_wide$Plantar_scales <- factor(disc_wide$Plantar_scales, levels = c("Smooth", "Pyramidal", "Spiny"), labels = c(1:length(setdiff(unique(disc_wide$Plantar_scales), NA))))
#disc_wide$Nasals <- factor(disc_wide$Nasals, levels = "contact", labels = c(1:length(setdiff(unique(disc_wide$Nasals), NA))))
#disc_wide$Prefrontal <- factor(disc_wide$Prefrontal, levels = c("separated", "point", "narrow", "contact"), labels = c(1:length(setdiff(unique(disc_wide$Prefrontal), NA))))

## Scale ordered discrete characters by maximum values:
rownames(disc_wide) <- disc_wide$REGNO ## Row names.
disc_wide <- disc_wide[2:ncol(disc_wide)]
for (var in names(disc_wide)) {
  disc_wide[[var]] <- as.numeric(disc_wide[[var]])
  disc_wide[[var]] <- disc_wide[[var]]/max(disc_wide[[var]], na.rm = TRUE) 
}
disc_wide$REGNO <- rownames(disc_wide) ## Add ID names back.

## Continuous variable and IDs:
cont_df <- morph_df[c("REGNO", "No_midbody_rows", "No_lamellae", "Snout_vent_length")] ## "SupLab" not variable.
rownames(cont_df) <- cont_df$REGNO ## Row names.
cont_df <- cont_df[2:ncol(cont_df)]

## Scale continuous variables by maximum values:
for (var in names(cont_df)) {
  cont_df[[var]] <- cont_df[[var]]/max(cont_df[[var]], na.rm = TRUE) 
}

## ID names back:
cont_df$REGNO <- rownames(cont_df)

## Combine discrete and continuous variables:
joint_df <- merge(cont_df, disc_wide, by = "REGNO")
NMDS_df <- joint_df[2:length(names(joint_df))] ## All characters, removing sample IDs.

## Replace missing data (i.e., NA) by mean value of each character:
for (char in names(NMDS_df)) { 
  NMDS_df[[char]] <- as.numeric(NMDS_df[[char]]) ## Convert to numeric.
  NMDS_df[[char]][is.na(NMDS_df[[char]])] <- mean(NMDS_df[[char]], na.rm = TRUE)
}

## Multidimensional scaling:
NMDS_out <- metaMDS(comm = NMDS_df, distance = 'euc', k = 2)

## Scores:
scores <- data.frame(REGNO = joint_df$REGNO,
                     Subspecies = morph_df$Taxon,
                     MDS1 = NMDS_out$points[, 1], 
                     MDS2 = NMDS_out$points[, 2])

## Combine with raw morph data:
NMDS_df$REGNO <- joint_df$REGNO
scores <- merge(NMDS_df, scores, by = "REGNO")

## Plot:
plot_biplot <- function(xvar, yvar) {

  ## Points:
  plot_NMDS <- ggplot(data = scores, aes(x = .data[[xvar]], y = .data[[yvar]])) +
  
    ## Plotting ellipses for each group:
    geom_mark_ellipse(data = filter(scores, Subspecies == "acripes"), color = palette[1], expand = unit(3, "mm"), size = 0.3, alpha = 0.25) +
    geom_mark_ellipse(data = filter(scores, Subspecies == "calx"), color = palette[2], expand = unit(3, "mm"), size = 0.3, alpha = 0.25) +
    geom_mark_ellipse(data = filter(scores, Subspecies == "ocellifer"), color = palette[3], expand = unit(3, "mm"), size = 0.3, alpha = 0.25) +
    geom_mark_ellipse(data = filter(scores, Subspecies == "pantherinus"), color = palette[4], expand = unit(3, "mm"), size = 0.3, alpha = 0.25) +
    
    ## Other params:
    #geom_point(aes(fill = Subspecies, shape = Subspecies), color = "black", size = 3, alpha = 0.75) +
    geom_jitter(aes(fill = Subspecies, shape = Subspecies), color = "black", size = 3, alpha = 0.75, width = 0.025, height = 0.025) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = shapes) +
    theme_bw() +  
    theme(#legend.position = c(0.1, 0.9),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 1, color = "gray20"),
          axis.title = element_text(size = 14, color = "gray20"),
          axis.text = element_text(size = 12, color = "gray20"),
          legend.title = element_text(size = 14, color = "gray20"),
          legend.text = element_text(size = 12, face = "italic", color = "gray20"),
          legend.key = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80"))
  
  ## Return:
  return(plot_NMDS)
  
} ## End of function.

## Run:
plot_biplot(xvar = "MDS1", yvar = "MDS2")
#plot_biplot(xvar = "SVL", yvar = "MBSR")
#plot_biplot(xvar = "Digital_lamellae", yvar = "Plantar_scales")

## PART 4: Combine plots ----

## NMDS plot:
plot_NMDS <- plot_biplot(xvar = "MDS1", yvar = "MDS2")
#ggsave(filename = paste0("figures/Fig_5_morphNMDS_2021-09-23.jpg"), plot = plot_NMDS, width = 15, height = 15, dpi = 500, units = "cm")
#ggsave(filename = paste0("figures/Fig_5_morphNMDS_2021-09-23.pdf"), plot = plot_NMDS, width = 24, height = 18, dpi = 500, units = "cm")

# ## Combine plots:
# plot_all <- ( plot_PCA | plot_NMDS ) + plot_annotation(tag_levels = "A") &
#   theme(plot.tag.position = c(0.01, 0.98), plot.tag = element_text(size = 22))
# plot_all
# ggsave(filename = paste0("figures/Fig_4_genPCA_morphNMDS_2021-08-09.jpg"), plot = plot_all, width = 25, height = 25, dpi = 500, units = "cm")

## PART 5: Combine clade info and morph data.

## Morph data:
char_df <- read.csv(file = "morphology/pantherinus_morphological_data_organized_ssps_added.csv", header = TRUE)
char_df <- char_df[c("REGNO", "Taxon", "Sex",
                     "No_midbody_rows", "No_lamellae", "Snout_vent_length", "SupLab",
                     "Dorsal_pattern", "Leg_pattern", "Digital_lamellae", "Plantar_scales", "Prefrontal")]
names(char_df)[names(char_df) == "REGNO"] <- "SAMPLE_ID"
names(char_df)[names(char_df) == "Taxon"] <- "SUBSPECIES"
char_df$Dorsal_pattern <- gsub(char_df$Dorsal_pattern, pattern = "Vert.", replacement = "Vert")
char_df$Dorsal_pattern <- gsub(char_df$Dorsal_pattern, pattern = " ", replacement = "_")
char_df$Digital_lamellae <- gsub(char_df$Digital_lamellae, pattern = " ", replacement = "_")
char_df$Prefrontal <- gsub(char_df$Prefrontal, pattern = "narrowly ", replacement = "")
char_df$Prefrontal <- gsub(char_df$Prefrontal, pattern = "narrow", replacement = "separated")
char_df$Sex <- gsub(char_df$Sex, pattern = "0", replacement = "m")
char_df$Sex <- gsub(char_df$Sex, pattern = "1", replacement = "f")
char_df$Sex <- gsub(char_df$Sex, pattern = "SA", replacement = "j")
char_df[is.na(char_df)] <- ""

## Clade info:
mit <- read.csv("data/assignments_mit_clades.csv", header = TRUE)
ddRAD <- read.csv("data/assignments_ddRAD_clades.csv", header = TRUE)
clades <- merge(mit, ddRAD, by = "SAMPLE_ID", all.x = TRUE, all.y = TRUE)
names(clades) <- c("SAMPLE_ID", "clade_mit", "clade_ddRAD")

## Subspecies info for sequenced specimens:
ssp <- read.csv("Supplementary_Material/Table_S1_genetic_samples.csv", header = TRUE)
ssp <- ssp[ssp$SPECIES == "pantherinus" & !is.na(ssp$SPECIES), ]
ssp <- arrange(ssp, "SAMPLE_ID")

## Merge and names:
clades <- merge(clades, ssp, by = "SAMPLE_ID", all.x = TRUE, all.y = TRUE)
clades <- clades[c("SAMPLE_ID", "clade_mit", "clade_ddRAD")]
clades$SAMPLE_ID <- gsub(clades$SAMPLE_ID, pattern = "_Ct_pant", replacement = "")
clades$SAMPLE_ID <- gsub(clades$SAMPLE_ID, pattern = "SAMR", replacement = "SAMAR")

## Merge with morph:
all_df <- merge(char_df, clades, by = "SAMPLE_ID", all.y = TRUE, all.x = TRUE)

## How many sequenced with morphology?
nrow(all_df[!is.na(all_df$No_midbody_rows), ])

## Save:
#write.csv(df, file = "data/matrix_morph-gen.csv", row.names = FALSE)

## PART 6: Summaries.

## Reload after minor edits:
#df <- read.csv("data/matrix_morph-gen.csv", header = TRUE)

## Discrete characters:
dis_df <- char_df[c(2, 8:12)]
molten <- melt(dis_df, id.vars = 1)
molten$var <- paste0(molten$variable, "_", molten$value)
dis_df <- dcast(molten, SUBSPECIES ~ var, length)
dis_df <- dis_df[setdiff(names(dis_df), c("Dorsal_pattern_", "Leg_pattern_", "Nasals_", "Plantar_scales_", "Prefrontal_"))]
dis_df$N <- table(char_df$SUBSPECIES)

## Calculate fractions of traits:
for (column in names(dis_df)[-1]) {
  dis_df[[column]] <- round(dis_df[[column]]/dis_df$N, 2)
}

## Columns:
dis_df <- dis_df[setdiff(names(dis_df), "N")]

## Continuous characters:
con_df <- char_df[c(2, 4:7)]
con_df[[2]] <- as.numeric(con_df[[2]])
con_df[[3]] <- as.numeric(con_df[[3]])
con_df[[4]] <- as.numeric(con_df[[4]])
con_df[[5]] <- as.numeric(con_df[[5]])

## Morph summaries by subspecies:
con_df <- con_df %>% group_by(SUBSPECIES) %>% summarise(N = length(SUBSPECIES),
                                                        mean_SVL = round(mean(Snout_vent_length, na.rm = TRUE), 1),
                                                        min_SVL = min(Snout_vent_length, na.rm = TRUE),
                                                        max_SVL = max(Snout_vent_length, na.rm = TRUE),
                                                        mean_nmr = round(mean(No_midbody_rows, na.rm = TRUE), 1),
                                                        min_nmr = min(No_midbody_rows, na.rm = TRUE),
                                                        max_nmr = max(No_midbody_rows, na.rm = TRUE),
                                                        mean_nla = round(mean(No_lamellae, na.rm = TRUE), 1),
                                                        min_nla = min(No_lamellae, na.rm = TRUE),
                                                        max_nla = max(No_lamellae, na.rm = TRUE))

## Means and ranges:
con_df$SVL <- paste0(con_df$mean_SVL, " (", con_df$min_SVL, "-", con_df$max_SVL, ")")
con_df$"No. midbody rows" <- paste0(con_df$mean_nmr, " (", con_df$min_nmr, "-", con_df$max_nmr, ")")
con_df$"No. digital lamellae" <- paste0(con_df$mean_nla, " (", con_df$min_nla, "-", con_df$max_nla, ")")
          
## Select:
con_df <- con_df[c(1, 2, 12, 13, 14)]
                                      
## Combine characters:
sum_df <- merge(con_df, dis_df, by = "SUBSPECIES")

## Save:
View(sum_df)
write.csv(sum_df, file = "data/morphology_summary_by_subspecies.csv", row.names = FALSE)

# End of script.
