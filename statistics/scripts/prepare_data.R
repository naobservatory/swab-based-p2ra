library(readr)
library(lubridate)
library(dplyr)
library(tidyr)
library(here)

# Helper script to load and preprocess data
# Sourced by by fit_model.R and quarto notebooks.

# Load data

tables_dir <- here("tables")

swab_metadata <- read_tsv(here(tables_dir, "swab-sample-metadata.tsv")) %>%
  mutate(date = ymd(paste0("20", date)))

# Use the per-day swab data, which disregards treatment
swab_reads <- read_tsv(here(tables_dir, "swabs-ra-summary.tsv")) %>%
  mutate(date = ymd(paste0("20", date))) %>%
  # Redundant with the metadata table
  select(-all_reads)

wastewater_metadata <- read_tsv(here(tables_dir, "wastewater-sample-metadata.tsv")) %>%
  mutate(date = ymd(paste0("20", date)))

wastewater_reads <- read_tsv(here(tables_dir, "ww-ra-summary.tsv")) %>%
  mutate(date = ymd(paste0("20", date))) %>%
  # Redundant with the metadata table
  select(-all_reads)

# Aggregate to species level and fill in zero counts for species unobserved in
# each sample
prepare_data <- function(all_species, metadata, read_data){
  all_species_and_samples <- crossing(all_species, metadata)
  # The taxonomic assignments go to below species level, but we want to analyze
  # at the species-level assignments
  read_data_species <- summarise(
    read_data,
    .by = c(date, location, species, group),
    viral_reads = sum(dedup_hv)
  )
  left_join(
    all_species_and_samples,
    read_data_species,
    by = join_by(date, location, species, group),
  ) %>%
    arrange(group, species) %>%
    mutate(viral_reads = replace_na(viral_reads, 0))
}

all_species <- bind_rows(swab_reads, wastewater_reads) %>%
  select(species, group) %>%
  distinct %>%
  arrange(group, species)

swab_reads_complete <- prepare_data(
  all_species,
  swab_metadata,
  swab_reads
  )
wastewater_reads_complete <- prepare_data(
  all_species,
  wastewater_metadata,
  wastewater_reads
  )
