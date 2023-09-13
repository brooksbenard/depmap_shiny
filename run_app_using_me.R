# runn_app_using_me.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 06/22/2023
# updated 09/13/2023
# Description: This script loads required packages, rds files, and helper funcitons required to run the DepMap Shiny app that performes mutation and cancer specific enrichment analyses based on a gene of interest's effect distribution.

# Set the working directory to the folder containing your Shiny app's files
setwd("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny")

# load packages ----
packages <-
  c(
    "scales",
    "ggplot2",
    "readxl",
    "reshape2",
    "plyr",
    "dplyr",
    "data.table",
    "stringr",
    "ggpubr",
    "janitor",
    "tidyverse",
    "magrittr",
    "cowplot",
    "ggExtra",
    "ggrepel",
    "tidyverse",
    "viridis",
    "cowplot",
    "RColorBrewer",
    "grid",
    "ggpubr",
    "gridExtra",
    "fgsea",
    "qusage",
    "shiny",
    "googledrive"
  )

# fgsea and qusage need to be downloaded through BiocManager
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("fgsea")
BiocManager::install("qusage")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# fastmap version needs to be specific for newer versions of R Studio
install.packages("fastmap", version = "1.1.1")

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# set global option and define "not in" function
`%ni%` <- Negate(`%in%`)
options(scipen = 999)

# load data ----
effect_scores <- 
  read_rds("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/effect_scores.rds")

mutations <-
  read_rds("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/mutations.rds")

# source functions ----
source("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/helpers.R")

# Run the Shiny app
runApp()
