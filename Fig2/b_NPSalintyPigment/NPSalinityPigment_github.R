# Plot scatterpie charts of phytoplankton community composition, derived from phytoplankton marker pigment ratios.
# Data: NASA S-MODE IOP-2 Bottle Data, CalCOFI bottle data (1954-2024)
# Links:
# https://podaac.jpl.nasa.gov/dataset/SMODE_L2_SHIPBOARD_BOTTLES_V1
# https://calcofi.org/data/oceanographic-data/bottle-database/
# Author: Shailja Gangrade
# Last updated: 9 July 2025

# Notes:

# The package 'PieGlyph' solves the problem of variable axes proportions;
# the package 'scatterpie' is only functional with equal coordinate axes as
# pies become distorted when changing plot dimensions. 

# Below, we use Tot_Chl_a (measured pigment in mg/m3) for pie radius size.
# Metadata on PO.DAAC lists Tot_Chl_a = monovinyl + divinyl chlorophyll a +
# allomers and epimers + chlorophyllide a.
# The data in chl_ug_L is simply chlorophyll concentration (ug/L)

# Ultimately pie chart sizes are not scaled well to the legend.
# So a manual post-export scaling of pie markers was applied using
# Adobe Illustrator tools. Contact author for more details.

################################################################################

## Load libraries
library(ggplot2)
library(scatterpie)
library(dplyr)
library(tidyr)
library(PieGlyph)

## Load data
# data_dir <- '/path/to/data/'
# fig_dir <- '/path/to/figures/'

btldata <- read.csv('IOP2-bottledata-pigmentratios-updatedMay2026.csv')
calcofi.data <- read.csv('CalCOFI_nutrient_data_updatedMay2026.csv')
calcofi.data.surf <- read.csv('CalCOFI_nutrient_data_surface_updatedMay2026.csv')

## Add column for N:P ratio
btldata$N.P.ratio <- btldata$nitrate_umol_L/btldata$phosphate_umol_L

## Replace NaN with 0 for pigment ratio data
btldata <- btldata %>%
  mutate(Perid.Chla = replace_na(Perid.Chla, 0),
         Allo.Chla = replace_na(Allo.Chla, 0),
         But.fuco.Chla = replace_na(But.fuco.Chla, 0),
         Fuco.Chla = replace_na(Fuco.Chla, 0),
         Hex.fuco.Chla = replace_na(Hex.fuco.Chla, 0),
         Zea.Chla = replace_na(Zea.Chla, 0),
         MVChlb.Chla = replace_na(MVChlb.Chla, 0))

## Subset for only near-surface data
btldata.surf <- btldata[btldata$depth<5,]

## Set up list of the selected pigment ratio variables and associated taxa names
pigment_ratio_vars <- c('Perid.Chla','Allo.Chla','But.fuco.Chla','Fuco.Chla',
                        'Hex.fuco.Chla','Zea.Chla','MVChlb.Chla')

pigment_vars <- c('Perid','Allo','ButFuco','Fuco',
                  'HexFuco','Zea','MVChlb')

taxa_names <- c('Dinoflagellate','Cryptophyte','Pelagophytes/Dictyochophytes',
                'Diatom','Haptophyte','Cyanobacteria','Chlorophyte')

## Set up color palette; derived'Paul Tol' color scheme

clrs <- c('#332288','#88CCEE','#44AA99','#117733','#DDCC77','#CC6677','#AA4499') # Manual selection from ptol colors.

################################################################################

## Plot pigment  composition pie charts in N:P vs. salinity space, with
## size of pie proportional to Total Chl-a concentration.

pigment_comp.surf <- ggplot(data=btldata.surf,aes(x=sal_abs, y=N.P.ratio)) +
  geom_point(data=calcofi.data.surf, aes(x=SA,y=NO3/PO4),color="grey",alpha=1, shape=17) +
  geom_hline(yintercept=16, linetype='dashed', col = 'black') +
  geom_pie_glyph(data=btldata.surf, aes(radius=Tot_Chl_a), slices=pigment_ratio_vars, alpha=0.6) +
  labs(x = "Salinity", y = "N:P", fill = "Pigment", radius = bquote('Chl-a (mg'~m^-3 ~')')) +
  scale_x_continuous(limits=c(32.5,34),breaks=seq(32.5,34,by=0.5)) +
  scale_y_continuous(limits=c(0,20),breaks=seq(0,20,by=5)) +
  scale_fill_manual(values=clrs,labels=pigment_vars)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))

pigment_comp.surf

################################################################################
