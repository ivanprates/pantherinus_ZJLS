###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To organize and make maps of morphological characters scored in museum subspecies of Ctenotus pantherinus.

## PART 1: Getting ready ----

## Packages:
library(patchwork)
library(raster)
library(reshape2)
library(rnaturalearth)
library(sf)
library(tidyverse)
library(scatterpie)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_pantherinus/")
detach("package:here", unload = TRUE)
library(here)

## Color palette:
palette <- c("#7AA5FB", "gray90", "#FFB000", "#5D5D5D", "#BB259C")

## Australia map:
AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113.25, -10.6, 153.7, -39.2), nrow = 2))) ## Extent.

## Base map:
base_map <- ggplot() +
  
  ## Plot base map:
  geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray40", size = 0.3) +
  
  ## Setting a bunch of parameters:
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank())     

## Check:
base_map

## Read shape file with taxon distribution polygons:
shape_panth <- read_sf(dsn = "range_shape", layer = "Ctenotus_pantherinus_range")

## Map of pantherinus range:
panth_map <- base_map + 
  geom_sf(data = shape_panth, fill = "blue", color = "transparent", size = 0.2, alpha = 0.5) +
  geom_sf(data = st_as_sf(AUS_map), fill = "transparent", color = "gray40", size = 0.4)
panth_map  

## Save:
#ggsave(filename = paste0("figures/pantherinus_range_map.pdf"), plot = panth_map, width = 15, height = 15, dpi = 500, units = "cm")

## PART 2: Character data ----

## Load and combine data:
WAM_df <- read.csv(file = "morphology/panther_scores_WAM_2021-09-17.csv", header = TRUE)
WAM_df$Rego <- paste0("WAMR_", WAM_df$Rego)
UMMZ_df <- read.csv(file = "morphology/panther_scores_UMMZ_2021-08-31.csv", header = TRUE)
UMMZ_df$Rego <- paste0(UMMZ_df$Museum, "_", UMMZ_df$Rego)
UMMZ_df <- UMMZ_df[setdiff(names(UMMZ_df), "Museum")]
char_df <- rbind(WAM_df, UMMZ_df)

## Size of juveniles (i.e., subadults):
char_df$SVL[char_df$Sex == "j"] <- NA
char_df$SVL[char_df$Sex == "J"] <- NA
char_df$SVL[char_df$Sex == "SA"] <- NA
char_df$SVL[char_df$SVL < 75] <- NA ## Examinations suggest specimens this small are juveniles.

## Average for subdigital lamellae counts:
char_df <- char_df %>% rowwise() %>% mutate(No_lamellae = mean(c(X4TLam1, X4TLam2), na.rm = TRUE))

## A few changes:
names(char_df)[names(char_df) == "Rego"] <- "REGNO"

char_df$Soles <- as.character(char_df$Soles)
char_df$Soles[char_df$Soles == "Rounded"] <- "Smooth"
char_df$Soles[char_df$Soles == "Pyramid"] <- "Pyramidal"
char_df$Soles[char_df$Soles == "Pryamid"] <- "Pyramidal"
char_df$Soles[char_df$Soles == "Spiky"] <- "Spiny"
char_df$Soles[char_df$Soles == "Spikey"] <- "Spiny"
char_df$Soles[char_df$Soles == "very spiky"] <- "Spiny"
char_df$Soles <- as.factor(char_df$Soles)

char_df$Dorsum <- as.character(char_df$Dorsum)
char_df$Dorsum[char_df$Dorsum == "Normal"] <- "Typical"
char_df$Dorsum[char_df$Dorsum == "Heavy"] <- "Wide border"
char_df$Dorsum[char_df$Dorsum == "Hiatus"] <- "Vert. hiatus"
char_df$Dorsum[char_df$Dorsum == "Hiatus, lines"] <- "Vert. hiatus"
char_df$Dorsum[char_df$Dorsum == "Some longitudinal lines"] <- "Lines"
char_df$Dorsum[char_df$Dorsum == "Vertebral stripe, half"] <- "Vert. stripe"
char_df$Dorsum[char_df$Dorsum == "Longitudinal lines"] <- "Lines"
char_df$Dorsum[char_df$Dorsum == "Vertebral stripe"] <- "Vert. stripe"
char_df$Dorsum <- as.factor(char_df$Dorsum)

char_df$Lamellae <- as.character(char_df$Lamellae)
char_df$Lamellae[char_df$Lamellae == "Broad"] <- "Single broad"
char_df$Lamellae[char_df$Lamellae == "Slighty broader"] <- "Single broad"
char_df$Lamellae[char_df$Lamellae == "Single"] <- "Single fine"
char_df$Lamellae[char_df$Lamellae == "Single mucron"] <- "Single fine"
char_df$Lamellae[char_df$Lamellae == "Single fine-ish"] <- "Single fine"
char_df$Lamellae[char_df$Lamellae == "Triple"] <- "Triple fine"
char_df$Lamellae[char_df$Lamellae == "Triple mucron"] <- "Triple fine"
char_df$Lamellae <- as.factor(char_df$Lamellae)

char_df$Legs <- as.character(char_df$Legs)
char_df$Legs[char_df$Legs == "lines"] <- "Lines"
char_df$Legs <- as.factor(char_df$Legs)

## Character names:
names(char_df)[names(char_df) == "Dorsum"] <- "Dorsal pattern"
names(char_df)[names(char_df) == "Lamellae"] <- "Digital lamellae"
names(char_df)[names(char_df) == "Soles"] <- "Plantar scales"
names(char_df)[names(char_df) == "MBSR"] <- "No. midbody rows"
names(char_df)[names(char_df) == "No_lamellae"] <- "No. lamellae"
names(char_df)[names(char_df) == "SVL"] <- "Snout-vent length"
names(char_df)[names(char_df) == "Legs"] <- "Leg pattern"

## Sample info, WAM:
WAMR_info <- read.csv("sample_information/WAM_Ct_pantherinus_Feb_2021.csv", header = TRUE)
WAMR_info$SITE <- as.character(WAMR_info$SITE)
WAMR_info <- WAMR_info[c("REGNO", "LATDEC", "LONGDEC", "SITE")]
WAMR_info$REGNO <- paste0("WAMR_", WAMR_info$REGNO)

## Sample info, all else:
UMMZ_info <- read.csv("sample_information/sample_information_pantherinus_ddRAD_2021-08-06.csv", header = TRUE)
UMMZ_info <- UMMZ_info[c("VOUCHER", "LAT", "LON", "LOCATION")]
names(UMMZ_info) <- c("REGNO", "LATDEC", "LONGDEC", "SITE")
UMMZ_info$REGNO <- gsub(UMMZ_info$REGNO, pattern = "SAMR", replacement = "SAMAR")
UMMZ_info <- UMMZ_info[UMMZ_info$REGNO %in% setdiff(UMMZ_info$REGNO, UMMZ_info$REGNO[grepl(UMMZ_info$REGNO, pattern = "WAMR")]) &
                       UMMZ_info$REGNO %in% setdiff(UMMZ_info$REGNO, c("no_voucher", NA))   , ]

## Merge:
sample_info <- rbind(WAMR_info, UMMZ_info)
char_df2 <- merge(char_df, sample_info, by = "REGNO", all.x = TRUE)

## Adding missing locality info:
char_df2$LATDEC[char_df2$REGNO == "UMMZ_242637"] <- -18.11428
char_df2$LONGDEC[char_df2$REGNO == "UMMZ_242637"] <- 123.5507
char_df2$SITE[char_df2$REGNO == "UMMZ_242637"] <- "Dampier Downs Station, Western Australia"

char_df2$LATDEC[char_df2$REGNO == "UMMZ_242640"] <- -19.23930
char_df2$LONGDEC[char_df2$REGNO == "UMMZ_242640"] <- 127.8264
char_df2$SITE[char_df2$REGNO == "UMMZ_242640"] <- "Carranya Station, Western Australia"

char_df2$LATDEC[char_df2$REGNO == "WAMR_177281"] <- -20.01261 
char_df2$LONGDEC[char_df2$REGNO == "WAMR_177281"] <- 120.86917
char_df2$SITE[char_df2$REGNO == "WAMR_177281"] <- "Sandfire airport"

char_df2$LATDEC[char_df2$REGNO %in% c("UMMZ_244307", "UMMZ_244310", "UMMZ_244312")] <-  -26.05799
char_df2$LONGDEC[char_df2$REGNO %in% c("UMMZ_244307", "UMMZ_244310", "UMMZ_244312")] <- 113.61562
char_df2$SITE[char_df2$REGNO %in% c("UMMZ_244307", "UMMZ_244310", "UMMZ_244312")] <- "PERON PENINSULA"

char_df2$LATDEC[char_df2$REGNO %in% c("UMMZ_244295", "UMMZ_244298")] <- -26.20080
char_df2$LONGDEC[char_df2$REGNO %in% c("UMMZ_244295", "UMMZ_244298")] <- 121.3036
char_df2$SITE[char_df2$REGNO %in% c("UMMZ_244295", "UMMZ_244298")] <- "LORNA GLEN STATION"

char_df2$LATDEC[char_df2$REGNO %in% c("UMMZ_210448", "UMMZ_210449")] <- -21.2255556000
char_df2$LONGDEC[char_df2$REGNO %in% c("UMMZ_210448", "UMMZ_210449")] <- 134.1300000000
char_df2$SITE[char_df2$REGNO %in% c("UMMZ_210448", "UMMZ_210449")] <- "49 KM S OF WAUCHOPE ON STUART HIGHWAY"

char_df2$SITE[char_df2$REGNO == "UMMZ_210450"] <- "27.1 KM S OF TI TREE ON STUART HIGHWAY"
char_df2$LATDEC[char_df2$REGNO == "UMMZ_210450"] <- -22.3608333
char_df2$LONGDEC[char_df2$REGNO == "UMMZ_210450"] <- 133.3877778

char_df2$SITE[char_df2$REGNO == "UMMZ_132495"] <- "6 MI S W OF WINNING TURNOFF ON NORTHWEST COASTAL HWY"
char_df2$LATDEC[char_df2$REGNO == "UMMZ_132495"] <- -23.218578
char_df2$LONGDEC[char_df2$REGNO == "UMMZ_132495"] <- 114.473225

## Rounding lat-longs:
char_df2$LATDEC <- round(char_df2$LATDEC, digits = 2)
char_df2$LONGDEC <- round(char_df2$LONGDEC, digits = 2)

## Combine sites very close to each other:
char_df2$SITE[char_df2$SITE == "AIRSTRIP SE CORNER OF BARROW ISLAND"] <- "BARROW ISLAND"
char_df2$SITE[char_df2$SITE == "SW CORNER OF BARROW ISLAND"] <- "BARROW ISLAND"
char_df2$SITE[char_df2$SITE == "12.5KM 309DEGREES (NEW) LISSADELL HS"] <- "AROUND NEW LISSADELL HOMESTEAD"
char_df2$SITE[char_df2$SITE == "10.8KM 253DEGREES (NEW) LISSADELL HS"] <- "AROUND NEW LISSADELL HOMESTEAD"
char_df2$SITE[char_df2$SITE == "15KM WSW NEW LISSADELL HOMESTEAD"] <- "AROUND NEW LISSADELL HOMESTEAD"
char_df2$SITE[char_df2$SITE == "11KM W NEW LISSADELL HOMESTEAD"] <- "AROUND NEW LISSADELL HOMESTEAD"
char_df2$SITE[char_df2$SITE == "27KM NNE TURKEY CREEK"] <- "AROUND NEW LISSADELL HOMESTEAD"
char_df2$SITE[char_df2$SITE == "8KM SW ONSLOW"] <- "ONSLOW AREA"
char_df2$SITE[char_df2$SITE == "Rangers' Head Quarters Karijini NP"] <- "Karijini NP"
char_df2$SITE[char_df2$SITE == "Murmur Kiwirrkurra IPA"] <- "Kiwirrkurra IPA"
char_df2$SITE[char_df2$SITE == "Lake Mackay Kiwirrkurra IPA"] <- "Kiwirrkurra IPA"
char_df2$SITE[char_df2$SITE == "8KM SE KUNUNURRA"] <- "AROUND KUNUNURRA"
char_df2$SITE[char_df2$SITE == "10KM W KUNUNURRA VICTORIA HIGHWAY"] <- "AROUND KUNUNURRA"
char_df2$SITE[char_df2$SITE == "25KM N KUNUNURRA 6KM ALONG CAVE SPRINGS ROAD"] <- "AROUND KUNUNURRA"
char_df2$SITE[char_df2$SITE == "16.8KM ENE BLACKSTONE"] <- "AROUND BLACKSTONE"
char_df2$SITE[char_df2$SITE == "18.4KM NE BLACKSTONE"] <- "AROUND BLACKSTONE"
char_df2$SITE[char_df2$SITE == "3.2KM N PUNGKULPIRRI WATERHOLE"] <- "AROUND PUNGKULPIRRI WATERHOLE"
char_df2$SITE[char_df2$SITE == "12.2KM N PUNGKULPIRRI WATERHOLE"] <- "AROUND PUNGKULPIRRI WATERHOLE"
char_df2$SITE[char_df2$SITE == "KUTJUNTARI ROCKHOLE"] <- "AROUND PUNGKULPIRRI WATERHOLE"
char_df2$SITE[char_df2$SITE == "REBECCA CREEK"] <- "AROUND PUNGKULPIRRI WATERHOLE"
char_df2$SITE[char_df2$SITE == "HAMELIN [DUMP]"] <- "Hamelin Station"
char_df2$SITE[char_df2$SITE == "7KM W HAMELIN"] <- "Hamelin Station"
char_df2$SITE[char_df2$SITE == "22KM NNW COBURN HS"] <- "Hamelin Station"
char_df2$SITE[char_df2$SITE == "23KM SSW BOORABBIN"] <- "NEAR BOORABBIN"
char_df2$SITE[char_df2$SITE == "16MI S KARALEE"] <- "NEAR BOORABBIN"
char_df2$SITE[char_df2$SITE == "15.5KM 086DEGREES TOOMEY HILLS"] <- "NEAR BOORABBIN"
char_df2$SITE[char_df2$SITE == "2KM N YARDIE CREEK THE WATERCOURSE"] <- "AROUND GIRALIA STATION"
char_df2$SITE[char_df2$SITE == "10.6 km N Exmouth"] <- "AROUND GIRALIA STATION"
char_df2$SITE[char_df2$SITE == "YARDIE CREEKTHE WATERCOURSE"] <- "AROUND GIRALIA STATION"
char_df2$SITE[char_df2$SITE == "GIRALIA STATION"] <- "AROUND GIRALIA STATION"
char_df2$SITE[char_df2$SITE == "49 KM S OF WAUCHOPE ON STUART HIGHWAY"] <- "AROUND BARROW CREEK"
char_df2$SITE[char_df2$SITE == "Barrow Creek"] <- "AROUND BARROW CREEK"
char_df2$SITE[char_df2$SITE == "MUGGON STATION"] <- "AROUND MUNGAWOLAGUDI CLAYPAN"
char_df2$SITE[char_df2$SITE == "6KM NE MUNGAWOLAGUDI CLAYPAN"] <- "AROUND MUNGAWOLAGUDI CLAYPAN"
char_df2$SITE[char_df2$SITE == "6.3KM 172DEGREES MOUNT PERCY"] <- "AROUND Carpenter Gap Windjana Gorge National Park"
char_df2$SITE[char_df2$SITE == "Carpenter Gap Windjana Gorge National Park"] <- "AROUND Carpenter Gap Windjana Gorge National Park"

## Save data:
char_df3 <- char_df2
names(char_df3) <- gsub(names(char_df3), pattern = " ", replacement = "_")
names(char_df3) <- gsub(names(char_df3), pattern = "-", replacement = "_")
names(char_df3) <- gsub(names(char_df3), pattern = "\\.", replacement = "")
write.csv(char_df3, file = "morphology/pantherinus_morphological_data_organized.csv", row.names = FALSE)

## PART 3: Character maps ----

## Data:
char_df3 <- read.csv(file = "morphology/pantherinus_morphological_data_organized.csv", header = TRUE)

## Maps:
plot_char_map <- function(char) {
  
  ## Testing:
  #char <- "Dorsal pattern" 
  #char <- "Digital lamellae"
  #char <- "Leg pattern"
  #char <- "Plantar scales"
  #char <- "Snout-vent length"
  
  ## Data:
  plot_df <- char_df2[c("REGNO", "LONGDEC", "LATDEC", "SITE", char)]
  names(plot_df)[names(plot_df) == char] <- "char"
  plot_df$char[plot_df$char == ""] <- NA
  plot_df <- plot_df[!is.na(plot_df$char), ]
  
  ## Quantitative characters:
  if(is.numeric(plot_df$char)) {
    
    ## Function: mean character state and variation by site:
    by_site <- function(site) {
    #site <- plot_df$SITE[1]
    site_df <- plot_df[plot_df$SITE == site, ]
    return(data.frame(SITE = site,
                      char = mean(site_df$char, na.rm = TRUE),
                      LONGDEC = mean(site_df$LONGDEC),
                      LATDEC = mean(site_df$LATDEC)))
    } ## End of function.
    
    ## Run:
    freqs_df <- map_df(unique(plot_df$SITE), by_site) 
    
  ## Qualitative characters:
  } else {
  
    ## Function: Character frequencies by site:
    get_char_site <- function(site) {
      #site <- plot_df$SITE[11]
      site_df <- plot_df[plot_df$SITE == site, ]
      site_tb <- table(site_df$char)/nrow(site_df)
      site_tb <- t(as.data.frame(as.matrix(site_tb)))
      site_df <- site_df[1, setdiff(names(site_df), c("REGNO", "char"))]
      site_df <- cbind(site_df, site_tb)
      site_df <- site_df[setdiff(names(site_df), "V1")]
      return(site_df)
    } ## End of function.
    
    ## Run:
    freqs_df <- map_df(unique(plot_df$SITE), get_char_site)
  
    ## Order of character states and palettes:
    ## palette_n <- c("#7AA5FB", "gray90", "#FFB000", "#5D5D5D", "#BB259C") ## b, lg, y, dg, p.
    if (char == "Dorsal pattern") {
      freqs_df  <- freqs_df[c(1:3, 5, 8, 7, 4, 6)]
      palette_n <- c("#FFB000", "#BB259C", "#5D5D5D", "#7AA5FB", "gray90") }
    if (char == "Leg pattern") {
      freqs_df  <- freqs_df[c(1:3, 6, 4, 5)]
      palette_n <- c("#FFB000", "#BB259C", "#7AA5FB") }
    if (char == "Digital lamellae") { 
      palette_n <- c("#FFB000", "#7AA5FB", "#BB259C") }
    if (char == "Plantar scales") {
      freqs_df  <- freqs_df[c(1:3, 5, 4, 6)]
      palette_n <- c("#FFB000", "#7AA5FB", "#BB259C") }
  
    } ## Close if statement.
  
  ## Plot map:
  char_map <- base_map +
    
    ## Setting a bunch of aesthetic parameters:
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, margin = margin(b = 5)))
    
  ## Plot data and change palette:
  ## Continuous characters:
  if(is.numeric(plot_df$char)) { 
    plot_f <- char_map + 
              ## Points:
              geom_point(data = freqs_df, aes(x = LONGDEC, y = LATDEC, fill = char), size = 3, color = "black", shape = 21, alpha = 0.9, stroke = 0.3) +
              scale_fill_gradient(name = char, low = "#ffffd9", high = "black") +
              theme(legend.position = "left", legend.direction = "vertical") ## Legend position.
  
  ## Discrete characters:  
  } else {
    plot_f <- char_map + 
               ## Pie charts:
               geom_scatterpie(data = freqs_df, aes(x = LONGDEC, y = LATDEC), 
                               pie_scale = 1.25, alpha = 0.85, color = "black",
                               cols = names(freqs_df)[4:ncol(freqs_df)]) +
              scale_fill_manual(values = palette_n, name = char) +
              theme(legend.position = "right", legend.direction = "vertical") ## Legend position.
  }
  
  ## Return:
  #plot_f
  return(plot_f)
  
} ## End of function.

## Discrete characters:
plot_A <- plot_char_map("Dorsal pattern") ; plot_A
plot_B <- plot_char_map("Digital lamellae") ; plot_B
plot_C <- plot_char_map("Plantar scales") ; plot_C
plot_G <- plot_char_map("Leg pattern") ; plot_G

## Continuous characters:
plot_D <- plot_char_map("Snout-vent length") ; plot_D
plot_E <- plot_char_map("No. midbody rows") ; plot_E
plot_F <- plot_char_map("No. lamellae") ; plot_F

## Horizontal:
#plot_char_h <- ( plot_A | plot_B | plot_C ) / ( plot_D | plot_E | plot_F )
#ggsave(filename = "figures/Fig_5_character_maps_2021-08-18_h.jpg", plot = plot_char_h, width = 25, height = 25, dpi = 300, units = "cm")
#ggsave(filename = "figures/Fig_5_character_maps_2021-08-18_h.pdf", plot = plot_char_h, width = 25, height = 25, dpi = 500, units = "cm")

## Vertical:
plot_char_v <- ( plot_spacer() | plot_A ) / ( plot_D | plot_G ) / ( plot_E | plot_C )
#ggsave(filename = "figures/Fig_3_character_maps_2021-09-17_mean_v.jpg", plot = plot_char_v, width = 25, height = 20, dpi = 500, units = "cm")
ggsave(filename = "figures/Fig_3U_2022-03-25.pdf", plot = plot_char_v, width = 25, height = 20, dpi = 500, units = "cm")

## A few more traits:
plot_char_y <- ( plot_F | plot_B )
#ggsave(filename = "figures/Fig_3_character_maps_2021-09-17_mean_v.jpg", plot = plot_char_v, width = 25, height = 20, dpi = 500, units = "cm")
ggsave(filename = "figures/Fig_3D_2022-03-25.pdf", plot = plot_char_y, width = 25, height = 7, dpi = 500, units = "cm")

## Plot locality labels:
char_t <- char_df2[c("LONGDEC", "LATDEC", "SITE")] 
char_t <- char_t[!duplicated(char_t$SITE), ]
plot_T <- plot_A + ggrepel::geom_text_repel(data = char_t, aes(x = LONGDEC, y = LATDEC, label = SITE), size = 2)
plot_T
#ggsave(filename = "figures/Fig_3_sites.pdf", plot = plot_T, width = 25, height = 20, dpi = 500, units = "cm")

## End of script.
