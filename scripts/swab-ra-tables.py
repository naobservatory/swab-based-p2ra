#!/usr/bin/env python3
# Standard library imports
import csv
import os
import json
from collections import defaultdict, Counter
from datetime import datetime
from taxonomy import load_taxonomy_names

from metadata_utils import first_level_mapping, second_level_mapping, is_date_in_range

# Constants and paths
validation_output_dir = "validation-output"
TABLE_DIR = "tables"

with open("deliveries.json") as f:
    deliveries = json.load(f)

target_deliveries = deliveries["swab-deliveries"]


# ---------------------------------------------------------------------
# Look‑ups
# ---------------------------------------------------------------------

taxid_names = load_taxonomy_names()

read_counts = {}
try:
    with open(os.path.join(TABLE_DIR, "n_reads_per_swab_sample.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            read_counts[row["sample"]] = int(row["reads"])
except FileNotFoundError:
    raise Exception("n_reads_per_swab_sample.tsv not found. Run fetch_readcounts.py first.")


date_loc_read_counts = Counter()
treatment_read_counts = Counter()
sample_treatment = {}

for delivery in target_deliveries:
    with open(
        os.path.join("..", "mgs-metadata", "deliveries", delivery, "metadata.tsv")
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            if row.get("demultiplexed") == "False":
                continue

            date = datetime.strptime(row["date"], "%Y-%m-%d")
            if not is_date_in_range(date):
                continue
            location = row["fine_location"]
            treatment = row["notes"]
            sample_treatment[sample] = treatment
            n_reads = read_counts.get(sample, 0)
            date_loc_read_counts[(date, location)]              += n_reads
            treatment_read_counts[(date, location, treatment)]  += n_reads

# ---------------------------------------------------------------------
# Counting HV reads (deduplicated and non-deduplicated)
# ---------------------------------------------------------------------

# Initialize data structures
samples = defaultdict(Counter)  # (date, location, pathogen) -> counts
treatment_samples = defaultdict(Counter) # (date, location, pathogen, treatment) -> counts

# Process classified reads
with open(os.path.join(validation_output_dir, "swabs-classified-all-reads.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        date = datetime.strptime(row["date"], "%Y-%m-%d")
        if not is_date_in_range(date):
            continue

        duplicate = row["is_duplicate"]
        location = row["loc"]
        pathogen = row["genome_name"]
        sample = row["sample"]
        treatment = sample_treatment[sample]
        if duplicate == "True":
            samples[(date, location, pathogen)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["non_dedup"] += 1
        else:
            samples[(date, location, pathogen)]["dedup"] += 1
            samples[(date, location, pathogen)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["dedup"] += 1

# Process non-validated reads
with open(os.path.join(validation_output_dir, "swabs-non-validated-reads.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        date = datetime.strptime(row["date"], "%Y-%m-%d")
        if not is_date_in_range(date):
            continue

        duplicate = row["is_duplicate"]
        location = row["loc"]
        taxid = row["taxid"]
        sample = row["sample_id"]
        pathogen = taxid_names[int(taxid)]
        treatment = sample_treatment[sample]

        if duplicate == "True":
            samples[(date, location, pathogen)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["non_dedup"] += 1
        else:
            samples[(date, location, pathogen)]["dedup"] += 1
            samples[(date, location, pathogen)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["non_dedup"] += 1
            treatment_samples[(date, location, pathogen, treatment)]["dedup"] += 1



# ---------------------------------------------------------------------
# Getting total reads
# ---------------------------------------------------------------------

# Get total reads for each grouping category
for (date, location, pathogen), counts in samples.items():
    n_reads = date_loc_read_counts[(date, location)]
    samples[(date, location, pathogen)]["all_reads"] = n_reads

for (date, location, pathogen, treatment), counts in treatment_samples.items():
    n_reads = treatment_read_counts[(date, location, treatment)]
    treatment_samples[(date, location, pathogen, treatment)]["all_reads"] = n_reads


# ---------------------------------------------------------------------
# Output results
# ---------------------------------------------------------------------

# Output results
with open(os.path.join(TABLE_DIR, "swabs-ra-summary.tsv"), "w") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow([
        "date",
        "location",
        "assignment",
        "species",
        "group",
        "non_dedup_hv",
        "dedup_hv",
        "all_reads"
    ])
    # Sort samples by date
    sorted_samples = sorted(samples.items(), key=lambda x: x[0][0])
    for (date, location, pathogen), data in sorted_samples:
        species = first_level_mapping(pathogen)
        group = second_level_mapping(species)
        date_str = date.strftime("%y%m%d")
        writer.writerow([
            date_str,
            location,
            pathogen,
            species,
            group,
            data.get("non_dedup", 0),
            data.get("dedup", 0),
            data["all_reads"],
        ])

with open(os.path.join(TABLE_DIR, "swabs-ra-per-treatment-summary.tsv"), "w") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow([
        "date",
        "location",
        "assignment",
        "species",
        "group",
        "treatment",
        "non_dedup",
        "dedup",
        "all_reads"
    ])
    # Sort samples by date
    sorted_samples = sorted(treatment_samples.items(), key=lambda x: x[0][0])

    for (date, location, pathogen, treatment), data in sorted_samples:
        species = first_level_mapping(pathogen)
        group = second_level_mapping(species)
        date_str = date.strftime("%y%m%d")
        writer.writerow([
            date_str,
            location,
            pathogen,
            species,
            group,
            treatment,
            data.get("non_dedup", 0),
            data.get("dedup", 0),data["all_reads"]])
