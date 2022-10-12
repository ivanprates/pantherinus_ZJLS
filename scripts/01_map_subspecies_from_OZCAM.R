###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To make maps of the subspecies of Ctenotus pantherinus based on museum records from the OZCAM database.

## PART 1: Getting ready ----

## Packages:
library(patchwork)
library(raster)
library(rnaturalearth)
library(sf)
library(tidyverse)

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_pantherinus/" ## Lab Linux.
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## Create folders:
dir.create(paste0(path, "data"))

## Shapes:
shapes <- c(22, 24, 21, 23)

## Color palette:
palette <- c("#ae5f23", "#f3e2c6", "#15b6ad", "#152868")
palette <- palette[c(1,3,2,4)]

## PART 2: Configure base map ----

## Australia map:
AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -45), nrow = 2))) ## Extent.

## Plot map:
base_map <- ggplot() +
  
  ## Plot base map:
  geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray20", size = 0.25) +
  
  ## Color of points:
  scale_fill_manual(values = palette) +
  
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

## PART 3: Map: WAMR subspecies ----

## WAMR sample information:
WAMR_info <- read.csv("sample_information/WAM_Ct_pantherinus_Feb_2021.csv", header = TRUE)
WAMR_info$SAMPLE_ID <- paste0("WAMR_", WAMR_info$REGNO, "_Ct_pant")
names(WAMR_info)[names(WAMR_info) == "REGNO"] <- "REGO"
names(WAMR_info)[names(WAMR_info) == "LATDEC"] <- "LAT"
names(WAMR_info)[names(WAMR_info) == "LONGDEC"] <- "LON"
WAMR_info$REGO <- factor(WAMR_info$REGO, levels = WAMR_info$REGO) 
WAMR_info <- WAMR_info[c("REGO", "LAT", "LON", "SUBSPECIES")]
WAMR_info$REGO <- paste0("WAM_", WAMR_info$REGO)

## ALA sample information:
ALA_info <- read.csv("sample_information/Ctenotus_pantherinus_subspecies_OZCAM_2021-08-04/Ctenotus_pantherinus_subspecies_OZCAM_2021-08-04.csv", header = TRUE)
ALA_info$Catalogue.Number <- as.character(ALA_info$Catalogue.Number)
for (ID in as.character(ALA_info$Catalogue.Number[ALA_info$Institution == "Western Australian Museum"])) {
  ALA_info$Catalogue.Number[ALA_info$Catalogue.Number == ID] <- paste0("WAM_", ID) }
ALA_info <- ALA_info[c("Catalogue.Number", "Decimal.latitude..WGS84.", "Decimal.longitude..WGS84.", "Subspecies")]
names(ALA_info) <- c("REGO", "LAT", "LON", "SUBSPECIES")
ALA_info$SUBSPECIES <- gsub(x = ALA_info$SUBSPECIES, pattern = "Ctenotus pantherinus ", replacement = "")
ALA_info$REGO <- gsub(x = ALA_info$REGO, pattern = "R\\.", replacement = "")
ALA_info$REGO <- gsub(x = ALA_info$REGO, pattern = "R", replacement = "")
ALA_info$REGO <- gsub(x = ALA_info$REGO, pattern = "\\..+", replacement = "")

## Combine:
subsp_info <- rbind(WAMR_info, ALA_info)

## Remove missing data:
subsp_info <- subsp_info[!is.na(subsp_info$LAT), ]
subsp_info <- subsp_info[!is.na(subsp_info$LON), ]

## Rounding coordinates to two decimal places (0.01 degrees equals a precision of 1,111 meters, i.e., roughly 1 km):
subsp_info$LAT <- round(subsp_info$LAT, digits = 1)
subsp_info$LON <- round(subsp_info$LON, digits = 1)

## Remove duplicates:
subsp_info <- subsp_info[!duplicated(subsp_info$REGO), ]
subsp_info <- arrange(subsp_info, by = REGO)

## How many specimens?
length(subsp_info$REGO)

## Unique latlons for map:
subsp_info <- subsp_info %>% group_by(LAT, LON) %>% sample_n(size = 1)

## Plot map:
map_subsp <-  base_map +
  
  ## Points:
  geom_point(data = subsp_info, aes(x = LON, y = LAT, fill = SUBSPECIES, shape = SUBSPECIES), 
             alpha = 0.5, color = "black", size = 3, stroke = 0.25, show.legend = FALSE) +
  
  ## Labels:
  labs(title = "Subspecies designation (WAM specimens)") +
  
  ## Symbols:
  scale_shape_manual(values = shapes) +
  
  ## Setting up facets:
  facet_wrap(~SUBSPECIES, nrow = 2) +
  
  ## Other stuff:
  theme(plot.margin = margin(t = 10, b = 0, l = 0, r = 15))

## Check:
map_subsp

## Saving plot (with cowplot):
#ggsave(filename = paste0("figures/subspecies_maps.jpg"), plot = map_subsp, width = 30, height = 8, dpi = 500, units = "cm")
ggsave(filename = paste0("figures/subspecies_maps_OZCAM.jpg"), plot = map_subsp, width = 15, height = 15, dpi = 500, units = "cm")
ggsave(filename = paste0("figures/subspecies_maps_OZCAM.pdf"), plot = map_subsp, width = 15, height = 15, dpi = 500, units = "cm")

## End of script.
