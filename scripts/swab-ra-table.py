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
validation_output_dir = "validation-output"
START_DATE = datetime(2025, 1, 7)
END_DATE = datetime(2025, 2, 25)

target_deliveries = [
    "NAO-ONT-20250120-Zephyr8",
    "NAO-ONT-20250127-Zephyr9",
    "NAO-ONT-20250213-Zephyr10",
    "NAO-ONT-20250213-Zephyr10-QC",
    "NAO-ONT-20250220-Zephyr11",
    "NAO-ONT-20250226-Zephyr10-QC2",
    "NAO-ONT-20250313-Zephyr12",
    "NAO-ONT-20250328-Zephyr12b"
]

# Functions
def is_date_in_range(date):
    return START_DATE <= date <= END_DATE

# Load names
taxid_names = load_taxonomy_names()

# Load counts
read_counts = {}
with open(os.path.join("outputs", "n_reads_per_swab_sample.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        read_counts[row["sample"]] = int(row["reads"])


# Create per-category counts and get treatment for each sample
date_loc_read_counts = defaultdict(int)
treatment_read_counts = defaultdict(int)
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
            location = row["fine_location"]
            treatment = row["notes"]
            date_loc_read_counts[(date, location)] += read_counts.get(sample, 0)
            sample_treatment[sample] = treatment
            treatment_read_counts[treatment] += read_counts.get(sample, 0)



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



# Get reads for each grouping category
for (date, location, pathogen), counts in samples.items():
    n_reads = date_loc_read_counts[(date, location)]
    samples[(date, location, pathogen)]["all_reads"] = n_reads

for (date, location, pathogen, treatment), counts in treatment_samples.items():
    
    n_reads = treatment_read_counts[treatment]
    treatment_samples[(date, location, pathogen, treatment)]["all_reads"] = n_reads

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

with open("swabs-ra-per-treatment-summary.tsv", "w") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["date", "location", "species", "pathogen", "treatment", "non_dedup", "dedup", "all_reads"])
    for (date, location, pathogen, treatment), data in sorted(treatment_samples.items()):
        species = first_level_mapping[pathogen]
        date_str = date.strftime("%y%m%d")
        writer.writerow([date_str, location, species, pathogen, treatment, data.get("non_dedup", 0), data.get("dedup", 0), data["all_reads"]])
