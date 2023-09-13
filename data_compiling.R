# Author : Brooks Benard, bbenard@stanford.edu
# Date: 06/09/2023
# updated 09/13/2023
# Description: This script downloads, compiles, and formats CRISPR gene effect,
# cell line mutations, gene expression, and metadata from the DepMap database

# load packages ----
packages <-
  c("dplyr",
    "tidyverse",
    "qusage",
    "readr",
    "fgsea",
    "janitor")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# set global options
options(scipen = 999)
options(timeout = 900)

# make directories ----
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/")
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data")
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/results")
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/results/plots")
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/results/data")

# Gene effect scores ----
# read in the gene effect scores from genome-wide CRISPR screening
effect_scores <-
  read.csv("https://figshare.com/ndownloader/files/40448555")

# remove everything but the gene name from the column headers
names(effect_scores) <-
  gsub("\\...*", "", names(effect_scores))

# make some duplicate gene names unique
colnames(effect_scores) <-
  make.unique(colnames(effect_scores))

# Cell line information ----
# all cell line mutations
mutations <-
  read.csv("https://figshare.com/ndownloader/files/40449653") |>
  dplyr::select(
    ProfileID,
    VariantType,
    VariantInfo,
    DNAChange,
    ProteinChange,
    HugoSymbol,
    CCLEDeleterious,
    CosmicHotspot
  ) |>
  subset(CCLEDeleterious == "True" | CosmicHotspot == "True")

# sample mapping identifiers
identifiers <-
  read.csv("https://figshare.com/ndownloader/files/40449635") |>
  dplyr::select(
    ModelID, ProfileID
  ) |>
  unique()

# sample metadata
metadata <-
  read.csv("https://figshare.com/ndownloader/files/40448834") |>
  dplyr::select(
    ModelID,
    StrippedCellLineName,
    DepmapModelType,
    OncotreeSubtype,
    OncotreePrimaryDisease,
    OncotreeLineage
  )

# combine the identifiers and metadata 
metadata <- 
  dplyr::left_join(identifiers, metadata,  by = "ModelID")

# add metadata to the mutations
mutations <-
  dplyr::left_join(mutations, metadata, by = "ProfileID")

# filter mutations
mutations <- mutations |>
  subset(ModelID %in% effect_scores$ModelID) # subset mutations to only those present in cell lines with CRISPR screening data
  
# copy number status ----
copy_number <-
  read.csv("https://figshare.com/ndownloader/files/40448840")

# add cell line identifier info to the cna file
names(copy_number)[names(copy_number) == 'X'] <- 'ModelID'
copy_number <- 
  left_join(copy_number, identifiers, by = "ModelID")

# Gene expression ----
rna <-
  read.csv("https://figshare.com/ndownloader/files/40449107")

# remove the white space between HUGO and ENSG from the header
colnames(rna) <-
  sub(" .*", "", colnames(rna))

# add column name for cell line name column
colnames(rna)[1] <- "ProfileID"

# wrangle the gene expression data into a usable form for DESeq2
rna <- rna |>
  setNames(sub("\\..*", "", names(rna))) |> # remove the ENSG identifier from the column names
  t() |>
  as.data.frame() |>
  row_to_names(row_number = 1) |>
  rownames_to_column() |>
  group_by(rowname) |>
  mutate(HugoID = paste(rowname, ".", 1:n(), sep = "")) |>
  ungroup() |>
  mutate(rowname = NULL) |>
  column_to_rownames(var = "HugoID") |>
  mutate_if(is_character, as.integer)

# GSEA pathway information ----
gmt <- "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/h.all.v2023.1.Hs.symbols.gmt"
pathways.hallmark <-
  gmtPathways(gmt)

# analysis functions require (1) gene effect data with cell line information labels, (2) cell line mutations with information labels, and (3) gene expression. Create and write out these data to rds files.
# 
# write out mutation + metadata file
write_rds(x = mutations, file = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/mutations.rds")

# write out copy number file
write_rds(x = copy_number, file = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/cna.rds")

# create the final effect score dataframe by adding the sample metadata
effect_scores <- 
  left_join(effect_scores, metadata, by = "ModelID")

# write out the gene effect scires + metadata file
write_rds(x = effect_scores, file = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/effect_scores.rds")

# write out the gene expression file
write_rds(x = rna, file = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/DepMap_Shiny/raw_data/effect_scores.rds")

