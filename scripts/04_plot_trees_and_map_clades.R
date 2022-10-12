###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To plot a phylogenetic tree based on the ddRAD data.

## PART 1: Getting ready ----

## Packages:
library(ape)
library(ggtree)
library(phytools)
library(patchwork)
library(raster)
library(rnaturalearth)
library(scico)
library(sf)
library(tidyverse)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_pantherinus/")
detach("package:here", unload = TRUE)
library(here)

## Color palette by clade:
#shapes <- c(22, 24, 21, 23)
shapes <- c(21, 22, 23, 24, 25)
palette <- c("black", "#b6e2dc", "#0f85a0", "#c5352a", "#ffae57") ## Black, light green, cyano, red, peach.

## Color palette by subspecies:
#shapes <- c(22, 24, 21, 23)
#palette <- c("#ae5f23", "#f3e2c6", "#15b6ad", "#152868") ; palette <- palette[c(1, 3, 2, 4)] 

## PART 2: Configure base map ----

## Australia map:
AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -45), nrow = 2))) ## Extent.

## Plot map:
base_map <- ggplot() +
  
  ## Plot base map:
  geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray20", size = 0.25) +
  
  ## Setting a bunch of parameters:
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(size = 1, color = "gray20"), 
    panel.background = element_rect(fill = "aliceblue"),
    
    ## Facet labels:
    strip.text = element_text(size = 14, face = "italic", color = "gray20", margin = margin(t = 0, r = 0, b = 5, l = 0)),
    strip.background = element_blank(),
    
    ## Titles:
    plot.title = element_blank(),
    #plot.title = element_text(size = 18, hjust = 0.5),
    
    ## Legend:
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10, face = "italic"),
    legend.key = element_rect(fill = "transparent"),
    #legend.position = c(0.40, 0.15),
    legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80"))  ## Box around legend.

## Check:
base_map

## PART 3: Mit tree ----

## Read tree:
tree <- read.tree(here("RAxML/pantherinus_n121_2021-06-14_codons.tre"))
tree <- drop.tip(tree, tip = c("WAMR_139286_Ct_rubi_Meentheena", "WAMR_146582_Ct_rubi_Burrup_Peninsula")) ## Drop outgroups.
tree <- drop.tip(tree, tip = c("ACRISP001_NA_Ct_pant", "WAMTR_1592_Ct_pant", "WAMTR_1708_Ct_pant")) ## Drop samples with no info to be found.
tree <- ladderize(tree) ## Ladderize.

## Simplified IDs:
tree$tip.label <- gsub(tree$tip.label, pattern = "(.+_.+_Ct_pant)_.+", replacement = "\\1")

## Changing IDs for consistency:
tree$tip.label <- gsub(tree$tip.label, pattern = "(.+)_NA_", replacement = "NA_\\1_")
sort(tree$tip.label)

## Define clades:
list_span <- list()
list_span$"mit clade 1" <- c("CUMV_14586_Ct_pant", "CUMV_14585_Ct_pant") 
list_span$"mit clade 2" <- c("WAMR_153984_Ct_pant", "WAMR_154151_Ct_pant")
list_span$"mit clade 3" <- c("SAMR_44860_Ct_pant", "WAMR_166402_Ct_pant")
#list_span$"mit clade 3" <- c("SAMR_44860_Ct_pant", "WAMR_137414_Ct_pant")
#list_span$"mit clade 4" <- c("SAMR_49197_Ct_pant", "WAMR_166402_Ct_pant")

## Function:
get_clades <- function(clade) {

  ## Node spanning defining samples:
  node <- findMRCA(tree = tree, type = "node", tips = list_span[[clade]])
  
  ## It's descendant nodes and tips:
  clade_f <- getDescendants(tree = tree, node = node)
  
  ## How many tips (i.e., samples) in clade?
  n_tips <- length(tree$tip.label)
  
  ## Descendant tips only:
  tips <- clade_f[clade_f <= n_tips]
  
  ## Corresponding tip labels:
  tip_labels <- data.frame(labels = tree$tip.label, n = 1:n_tips) ## Creating data frame to make things easier.
  tips_clade <- tip_labels[tip_labels$n %in% tips, ]$labels ## Extracting tip labels in focal clade.
  
  ## return:
  return(data.frame(SAMPLE_ID = tips_clade, clade = clade, stringsAsFactors = FALSE))

} ## End of function.

## Run:
clade_df <- purrr::map_df(names(list_span), get_clades)

## Save:
#write.csv(clade_df, file = "data/assignments_mit_clades.csv", row.names = FALSE)

## Add lat-lon info:
sample_info <- read.csv(file = "sample_information/sample_information_pantherinus_cytb_2021-08-06.csv", header = TRUE)

## Combine:
clade_ll <- left_join(clade_df, sample_info, by = "SAMPLE_ID")

## Unique latlons for map:
clade_u <- clade_ll %>% group_by(LOCATION) %>% sample_n(size = 1)

## Column name:
names(clade_u)[names(clade_u) == "SUBSPECIES"] <- "Subspecies"

## Plot map:
map_mit <-  base_map +
  
  ## Points:
  geom_point(data = clade_u, aes(x = LON, y = LAT, fill = clade, shape = clade), 
             color = "black", size = 2.75, alpha = 0.85, stroke = 0.35, show.legend = FALSE) +
  
  ## Colors:
  scale_fill_manual(values = palette[c(3:5)]) +
  
  ## Plot acripes again on top for better visibility among ocellifer samples in the west:
  #geom_point(data = filter(clade_u, Subspecies == "acripes"), aes(x = LON, y = LAT, fill = Subspecies, shape = Subspecies), 
  #           color = "black", size = 2.75, alpha = 0.70, stroke = 0.25, show.legend = FALSE) +
  
  ## Labels:
  #ggrepel::geom_text_repel(data = clade_u, aes(x = LON, y = LAT, label = SAMPLE_ID), size = 2) +
  
  ## Labels:
  theme(strip.text = element_text(face = "plain")) +
  scale_shape_manual(values = shapes) +
  labs(fill = "Clade", shape = "Clade") #+
  
  ## Setting up facets:
  #facet_wrap(~clade, ncol = 1)

## Check:
map_mit

## Create list and group samples by clade to color tree tips:
clade_c <- clade_ll[c("SAMPLE_ID", "clade")] ## Simplifying.
clade_list <- split(clade_c, f = clade_c$clade) ## Split dataframe into list.
for (clade in names(clade_list)) { clade_list[[clade]] <- clade_list[[clade]]$SAMPLE_ID } ## Keeping only IDs.
tree <- groupOTU(tree, clade_list) ## Group samples to plot tree.

## Alternatively, group by subspecies designation:
clade_s <- clade_ll[c("SAMPLE_ID", "SUBSPECIES")] ## Simplifying.
subsp_list <- split(clade_s, f = clade_s$SUBSPECIES) ## Split dataframe into list.
for (subsp in names(subsp_list)) { subsp_list[[subsp]] <- subsp_list[[subsp]]$SAMPLE_ID } ## Keeping only IDs.
#tree <- groupOTU(tree, subsp_list) ## Group samples to plot tree.

## Prepare tree to plot:
tree_plot <- tree
tree_plot$node.label <- as.numeric(tree_plot$node.label)
tree_plot$node.label[tree_plot$node.label < 70] <- "" ## Remove support for unsupported nodes.

## Editing some tip labels:
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "NA_", replacement = "") 
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "_Ct_pant", replacement = "")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "_", replacement = " ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "CCM", replacement = "CCM ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "DLR", replacement = "DLR ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "PANTSP", replacement = "PANTSP ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "PMO", replacement = "PMO ")

## Nodes to add bars:
node_list <- list_span
for (clade in names(node_list)) { 
  node_list[[clade]] <- gsub(node_list[[clade]], pattern = "_Ct_pant", replacement = "")
  node_list[[clade]] <- gsub(node_list[[clade]], pattern = "_", replacement = " ")
  node_list[[clade]] <- findMRCA(tree = tree_plot, type = "node", tips = node_list[[clade]])
}

## Plot tree:
tree_mit <- ggtree(tree_plot, size = 0.5, color = "gray20", ladderize = FALSE) + ## Size = branch line thickness.
  #aes(color = group) ## If coloring tips by group.
  
  ## Editing tree tips:
  geom_tiplab(color = "gray20", size = 1.8, offset = 0.004) +
  #geom_tiplab(aes(color = group), size = 1.75, offset = 0.004, fontface = 2) + ## Color tips by group.
  
  ## If using clades:
  geom_tippoint(aes(shape = group, fill = group), color = "gray20", size = 1.5, stroke = 0.35) +
  scale_color_manual(values = palette[c(3:5)]) +
  scale_fill_manual(values = palette[c(3:5)]) +
  scale_shape_manual(values = c(21:23)) +
  
  ## If using subspecies:
  #geom_tippoint(aes(shape = group, fill = group), color = "gray20", size = 2.25, stroke = 0.5) +
  #scale_color_manual(values = palette, name = "Subspecies") +
  #scale_fill_manual(values = palette, name = "Subspecies") +
  #scale_shape_manual(values = shapes, name = "Subspecies") +
  
  ## If using clades, add clade bars:    
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.005, offset = 0.04, node = node_list[[1]], label = names(node_list)[[1]], color = c(palette[3], "gray20")) + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.005, offset = 0.04, node = node_list[[2]], label = names(node_list)[[2]], color = c(palette[4], "gray20")) + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.005, offset = 0.04, node = node_list[[3]], label = names(node_list)[[3]], color = c(palette[5], "gray20")) + 
  
  ## Node support:
  #geom_nodelab(size = 1.5, col = "gray50", nudge_x = -0.005, nudge_y = -0.75) +
  
  ## Other edits:
  theme(legend.position = "none")
  # theme(legend.position = c(0.85, 0.12),
  #       legend.title = element_text(size = 14, color = "gray20"),
  #       legend.text = element_text(size = 12, face = "italic", color = "gray20"),
  #       legend.key = element_rect(fill = "transparent"),
  #       legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80"))

## Plot:
tree_mit <- tree_mit + xlim(0, 0.23) ## By clade.
#tree_mit ## By subspecies.

## PART 4: Combine mit plots ----

## Combine plots:
#plot_mit <- ( tree_mit | map_mit ) + plot_layout(widths = c(2.5, 1))

## Save:
#ggsave(plot = plot_mit, filename = "figures/Fig_2_pantherinus_mit.jpg", width = 25, height = 24, units = "cm", dpi = 500)
#ggsave(plot = plot_mit, filename = "figures/Fig_2_pantherinus_mit.pdf", width = 25, height = 24, units = "cm", dpi = 500)

## PART 5: ddRAD tree ----

## Read tree:
#tree <- read.tree(here("RAxML/RAxML_bipartitions_ddRAD_c85_2021-06-14.tre")) ## 1st submission, SNPs only.
tree <- read.tree(here("RAxML/RAxML_bipartitions_ddRAD_c85_complete_matrix_2021-10-12.tre")) ## 2nd submission, included invariant sites.
tree <- ladderize(tree)

## Drop missidentified sample and outgroups:
tree <- drop.tip(phy = tree, tip = "WAMR_139249_Ct_saxa")
tree <- drop.tip(phy = tree, tip = tree$tip.label[grep(x = tree$tip.label, pattern = "_Ct_pant", invert = TRUE)])

## Define clades:
list_span <- list()
list_span$"nuc clade 1" <- "NA_PMO179_Ct_pant"
list_span$"nuc clade 2" <- c("SAMR_65403_Ct_pant", "SAMR_55674_Ct_pant")
#list_span$"nuc clade 1" <- c("NA_PMO179_Ct_pant", "SAMR_55674_Ct_pant")
list_span$"nuc clade 3" <- c("CUMV_14576_Ct_pant", "WAMR_129863_Ct_pant") 
#list_span$"nuc clade 4" <- c("WAMR_125012_Ct_pant", "UMMZ_242643_Ct_pant") ## 1st submission (SNPs only)
#list_span$"nuc clade 5" <- c("WAMR_140704_Ct_pant", "CUMV_14545_Ct_pant") ## 1st submission (SNPs only)
list_span$"nuc clade 4" <- c("WAMR_140704_Ct_pant", "UMMZ_242643_Ct_pant")
list_span$"nuc clade 5" <- c("SAMR_45528_Ct_pant", "CUMV_14545_Ct_pant")

## Function:
get_clades <- function(clade) {
  
  ## Node spanning defining samples:
  node <- findMRCA(tree = tree, type = "node", tips = list_span[[clade]])
  
  ## It's descendant nodes and tips:
  clade_f <- getDescendants(tree = tree, node = node)
  
  ## How many tips (i.e., samples) in clade?
  n_tips <- length(tree$tip.label)
  
  ## Descendant tips only:
  tips <- clade_f[clade_f <= n_tips]
  
  ## Corresponding tip labels:
  tip_labels <- data.frame(labels = tree$tip.label, n = 1:n_tips) ## Creating data frame to make things easier.
  tips_clade <- tip_labels[tip_labels$n %in% tips, ]$labels ## Extracting tip labels in focal clade.
  
  ## return:
  return(data.frame(SAMPLE_ID = tips_clade, clade = clade, stringsAsFactors = FALSE))
  
} ## End of function.

## Run:
clade_df <- purrr::map_df(c("nuc clade 2", "nuc clade 3", "nuc clade 4", "nuc clade 5"), get_clades)
clade_df <- rbind(clade_df, c("NA_PMO179_Ct_pant", "nuc clade 1"))
clade_df <- arrange(clade_df, by = clade)

## Save:
#write.csv(clade_df, file = "data/assignments_ddRAD_clades.csv", row.names = FALSE)

## Add lat-lon info:
sample_info <- read.csv(file = "sample_information/sample_information_pantherinus_ddRAD_2021-08-06.csv", header = TRUE)
clade_ll <- merge(sample_info, clade_df, by = "SAMPLE_ID")

## Unique latlons for map:
clade_u <- clade_ll %>% group_by(LOCATION) %>% sample_n(size = 1)

## Column name:
names(clade_u)[names(clade_u) == "SUBSPECIES"] <- "Subspecies"

## Plot map:
map_ddRAD <-  base_map +
  
  ## Points:
  geom_point(data = clade_u, aes(x = LON, y = LAT, fill = clade, shape = clade), color = "black", size = 2.75, 
             alpha = 0.85, stroke = 0.35, show.legend = FALSE) +
  
  ## Colors:
  scale_fill_manual(values = palette) +
  
  ## Labels:
  #ggrepel::geom_text_repel(data = clade_u, aes(x = LON, y = LAT, label = SAMPLE_ID), size = 2) +
  
  ## Labels:
  theme(strip.text = element_text(face = "plain")) +
  #scale_shape_manual(values = shapes) + ## By subspecies.
  scale_shape_manual(values = rev(shapes[c(3, 2, 1, 4, 5)])) + ## By clades.
  labs(fill = "Subspecies", shape = "Subspecies")

  ## Setting up facets:
  #facet_wrap(~clade, ncol = 1)
  
## Check:
map_ddRAD

## Create list and group samples by clade to color tree tips:
clade_s <- clade_df[c("SAMPLE_ID", "clade")] ## Simplifying.
clade_list <- split(clade_s, f = clade_s$clade) ## Split dataframe into list.
for (clade in names(clade_list)) { clade_list[[clade]] <- clade_list[[clade]]$SAMPLE_ID } ## Keeping only IDs.
tree <- groupOTU(tree, clade_list) ## Group samples to plot tree.

## Alternatively, group by subspecies designation:
clade_ll <- clade_ll[c("SAMPLE_ID", "SUBSPECIES")]
subsp_list <- split(clade_ll, f = clade_ll$SUBSPECIES) ## Split dataframe into list.
for (subsp in names(subsp_list)) { subsp_list[[subsp]] <- subsp_list[[subsp]]$SAMPLE_ID } ## Keeping only IDs.
#tree <- groupOTU(tree, subsp_list) ## Group samples to plot tree.

## Prepare tree to plot:
tree_plot <- tree
tree_plot$node.label <- as.numeric(tree_plot$node.label)
tree_plot$node.label[tree_plot$node.label < 70] <- "" ## Remove support for unsupported nodes.

## Editing some tip labels:
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "_NA_Ct_pant", replacement = "") 
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "_Ct_pant", replacement = "")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "_", replacement = " ")

## A few more edits:
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "NA DLR", replacement = "DLR ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "NA CCM", replacement = "CCM ")
tree_plot$tip.label <- gsub(tree_plot$tip.label, pattern = "NA PMO", replacement = "PMO ")

## Nodes to add bars:
node_list <- list_span
for (clade in names(node_list)) { 
  node_list[[clade]] <- gsub(node_list[[clade]], pattern = "_Ct_pant", replacement = "")
  node_list[[clade]] <- gsub(node_list[[clade]], pattern = "_", replacement = " ")
  node_list[[clade]] <- gsub(node_list[[clade]], pattern = "NA PMO", replacement = "PMO ")
  node_list[[clade]] <- findMRCA(tree = tree_plot, type = "node", tips = node_list[[clade]])
}

## Plot tree:
tree_ddRAD <- ggtree(tree_plot, size = 0.5, ladderize = FALSE, color = "gray20") + ## Size = branch line thickness.
  #aes(color = group) ## If coloring tips by group.
  
  ## Editing tree tips:
  geom_tiplab(color = "gray20", size = 2.25, offset = 0.0005) +
  #geom_tiplab(aes(color = group), size = 2.25, offset = 0.004, fontface = 2) + ## Color tips by group.
  geom_tippoint(aes(shape = group, fill = group), color = "gray20", size = 2.5, stroke = 0.5) +
  
  ## By subspecies:
  #scale_color_manual(values = palette, name = "Subspecies") +
  #scale_fill_manual(values = palette, name = "Subspecies") +
  #scale_shape_manual(values = shapes, name = "Subspecies") +
  
  ## By clades:
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_shape_manual(values = rev(shapes[c(3, 2, 1, 4, 5)])) +
  
  ## Clade bars:
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.0006, offset = 0.004, node = node_list[[1]], label = names(node_list)[[1]], color = c(palette[2], "gray20")) + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.0006, offset = 0.004, node = node_list[[2]], label = names(node_list)[[2]], color = c(palette[3], "gray20")) + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.0006, offset = 0.004, node = node_list[[3]], label = names(node_list)[[3]], color = c(palette[4], "gray20")) + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.0006, offset = 0.004, node = node_list[[4]], label = names(node_list)[[4]], color = c(palette[5], "gray20")) + 
  
  ## Node support:
  #geom_nodelab(size = 1.5, col = "gray50", nudge_x = -0.005, nudge_y = -0.75) +
  
  ## Other edits:
  theme(legend.position = "none") ## By clade.
  #theme(#legend.position = c(0.8, 0.11),
  #      legend.title = element_text(size = 14, color = "gray20"),
  #      legend.text = element_text(size = 12, face = "italic", color = "gray20"),
  #      legend.key = element_rect(fill = "transparent"),
  #      legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80"))

## Plot limits:
tree_ddRAD <- tree_ddRAD + xlim(0, 0.022) ## By clade.
#tree_ddRAD ## By subspecies.

## PART 6: Combine trees and maps based on clades ----

## Combine:
plot_maps <- (map_mit / map_ddRAD )
plot_h <- ( tree_mit | tree_ddRAD | plot_maps) + plot_layout(widths = c(1, 1, 0.75)) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.009, 0.95), plot.tag = element_text(size = 22))
plot_h

## Save:
#ggsave(plot = plot_h, filename = "figures/Fig_3_pantherinus_clades_h.jpg", width = 35, height = 20, units = "cm", dpi = 500)
#ggsave(plot = plot_h, filename = "figures/Fig_3_pantherinus_clades_h.pdf", width = 35, height = 20, units = "cm", dpi = 500)

## Combine:
plot_mit <- ( tree_mit / map_mit ) + plot_layout(heights =c (2, 1))
plot_ddRAD <- ( tree_ddRAD / map_ddRAD ) + plot_layout(heights = c(2, 1))
plot_v <- ( plot_mit | plot_ddRAD ) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0.08, 0.97), plot.tag = element_text(size = 22))

## Save:
#ggsave(plot = plot_v, filename = "figures/Fig_4.jpg", width = 25, height = 28, units = "cm", dpi = 500)
#ggsave(plot = plot_v, filename = "figures/Fig_4.pdf", width = 25, height = 28, units = "cm", dpi = 500)

## PART 7: Combine trees based on subspecies ----

## Combine plots:
#plot_ssp <- ( tree_mit | tree_ddRAD )

## Save:
#ggsave(plot = plot_ssp, filename = "figures/Fig_5.jpg", width = 25, height = 25.5, units = "cm", dpi = 500)
#ggsave(plot = plot_ssp, filename = "figures/Fig_5.pdf", width = 25, height = 25.5, units = "cm", dpi = 500)

## End of script.
