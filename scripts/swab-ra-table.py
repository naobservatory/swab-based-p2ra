#!/usr/bin/env python3
# Standard library imports
import csv
import json
import os
import sys
from collections import defaultdict, Counter
from datetime import datetime
from taxonomy import load_taxonomy_names

# Constants and paths
dashboard_dir = os.path.expanduser("~/code/mgs-restricted/dashboard")
validation_output_dir = "validation-output"
START_DATE = datetime(2025, 1, 7)
END_DATE = datetime(2025, 2, 25)

# Functions
def is_date_in_range(date):
    return START_DATE <= date <= END_DATE

# Load data
taxid_names = load_taxonomy_names()

read_counts = {}
with open("n_reads_per_swab_sample.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        read_counts[row["sample"]] = int(row["reads"])


# Create virus clade groupings
first_level_mapping = {
    # Coronaviruses (seasonal)
    "Human coronavirus OC43": "HCoV-OC43",
    "Human coronavirus 229E": "HCoV-229E",
    "Human coronavirus HKU1": "HCoV-HKU1",
    "Human coronavirus NL63": "HCoV-NL63",

    # Coronaviruses (SARS-CoV-2)
    "Severe acute respiratory syndrome coronavirus 2": "SARS-CoV-2",

    # Influenza
    "H1N1": "Influenza A",

    # Rhinoviruses
    "Rhinovirus A34": "Rhinovirus A",
    "Rhinovirus A94": "Rhinovirus A",
    "Rhinovirus B37": "Rhinovirus B",
    "Rhinovirus C2": "Rhinovirus C",
    "Rhinovirus C36": "Rhinovirus C",
    "Rhinovirus C42": "Rhinovirus C",
    "Rhinovirus C56": "Rhinovirus C",

    # Mononegavirales
    "RSVA": "RSV-A",
    "RSVB": "RSV-B",
    "HPIV4": "HPIV4b"
}


# Initialize data structures
samples = defaultdict(Counter)  # (date, location, pathogen) -> counts
date_loc_sample_map = defaultdict(set)  # (date, location) -> set of samples

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

        date_loc_sample_map[(date, location)].add(sample)
        if duplicate == "True":
            samples[(date, location, pathogen)]["non_dedup"] += 1
        else:
            samples[(date, location, pathogen)]["dedup"] += 1
            samples[(date, location, pathogen)]["non_dedup"] += 1


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

        date_loc_sample_map[(date, location)].add(sample)
        if duplicate == "True":
            samples[(date, location, pathogen)]["non_dedup"] += 1
        else:
            samples[(date, location, pathogen)]["dedup"] += 1
            samples[(date, location, pathogen)]["non_dedup"] += 1

# Calculate total reads
for (date, location, pathogen), counts in samples.items():
    n_reads = 0
    for sample in date_loc_sample_map[(date, location)]:
        n_reads += read_counts[sample]
    samples[(date, location, pathogen)]["all_reads"] = n_reads

# Output results
with open("swabs-ra-summary.tsv", "w") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["date", "location", "species", "assignment", "non_dedup_hv", "dedup_hv", "all_reads"])
    # Sort samples by date
    sorted_samples = sorted(samples.items(), key=lambda x: x[0][0])
    for (date, location, pathogen), data in sorted_samples:
        species = first_level_mapping[pathogen]
        date_str = date.strftime("%y%m%d")
        writer.writerow([
            date_str,
            location,
            species,
            pathogen,
            data.get("non_dedup", 0),
            data.get("dedup", 0),
            data["all_reads"],
        ])