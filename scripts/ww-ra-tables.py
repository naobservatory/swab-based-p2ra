#!/usr/bin/env python3
# Standard library imports
import csv
import json
import os
from collections import defaultdict, Counter
from datetime import datetime
from taxonomy import load_taxonomy_names

from metadata_utils import first_level_mapping, second_level_mapping, is_date_in_range, pathogens_to_ignore, parse_count_table

# Constants and paths
validation_output_dir = "validation-output"
TABLE_DIR = "tables"


with open("deliveries.json") as f:
    deliveries = json.load(f)

target_deliveries = deliveries["ww-deliveries"]

# ---------------------------------------------------------------------
# Look‑ups
# ---------------------------------------------------------------------

taxid_names = load_taxonomy_names()

read_counts = parse_count_table("n_reads_per_ww_sample")

date_loc_read_counts = Counter()

# Reading metadata to get total read counts per date-location
for delivery in target_deliveries:
    with open(
        os.path.join("..", "mgs-metadata", "deliveries", delivery, "metadata.tsv")
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            fine_location = row["fine_location"]
            if fine_location not in ("DNI", "DSI"):
                continue

            date = datetime.strptime(row["date"], "%Y-%m-%d")
            if not is_date_in_range(date):
                continue

            date = datetime.strptime(row["date"], "%Y-%m-%d")
            date_loc_read_counts[(date, fine_location)] += read_counts.get(sample, 0)

# ---------------------------------------------------------------------
# Counting HV reads (deduplicated and non-deduplicated)
# ---------------------------------------------------------------------

# Initialize data structures
samples = defaultdict(Counter)  # (date, location, pathogen) -> counts

# Process classified reads
with open(os.path.join(validation_output_dir, "ww-classified-all-reads.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        date = datetime.strptime(row["date"], "%Y-%m-%d")
        if not is_date_in_range(date):
            continue


        duplicate = row["is_duplicate"]

        key = date, row["loc"], row["genome_name"]
        if duplicate == "True":
            samples[key]["dedup"] += 1
        else:
            samples[key]["non_dedup"] += 1
            samples[key]["dedup"] += 1

# Process non-validated reads
with open(os.path.join(validation_output_dir, "ww-non-validated-reads.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        date = datetime.strptime(row["date"], "%Y-%m-%d")
        if not is_date_in_range(date):
            continue


        taxid = row["taxid"]
        pathogen = taxid_names[int(taxid)]

        duplicate = row["is_duplicate"]

        key = date, row["loc"], pathogen
        if duplicate == "True":
            samples[key]["dedup"] += 1
        else:
            samples[key]["non_dedup"] += 1
            samples[key]["dedup"] += 1



# ---------------------------------------------------------------------
# Getting total reads
# ---------------------------------------------------------------------

for (date, location, pathogen), counts in samples.items():
    n_reads = date_loc_read_counts[(date, location)]
    samples[(date, location, pathogen)]["all_reads"] = n_reads

# ---------------------------------------------------------------------
# Output results
# ---------------------------------------------------------------------

with open(os.path.join(TABLE_DIR, "ww-ra-summary.tsv"), "w") as outf:
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
    sorted_samples = sorted(samples.items())
    for (date, location, pathogen), data in sorted_samples:
        if pathogen in pathogens_to_ignore():
            continue
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