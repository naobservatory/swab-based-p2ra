#! /usr/bin/env python3
import pandas as pd
import csv
import json
import os
from datetime import datetime
from collections import defaultdict

validation_output_dir = "validation-output"

# Create virus clade groupings
clade_mapping = {
    # Coronaviruses (seasonal)
    "Human coronavirus OC43": "Coronaviruses (seasonal)",
    "Human coronavirus 229E": "Coronaviruses (seasonal)",
    "Human coronavirus HKU1": "Coronaviruses (seasonal)",
    "Human coronavirus NL63": "Coronaviruses (seasonal)",

    # Coronaviruses (SARS-CoV-2)
    "Severe acute respiratory syndrome coronavirus 2": "Coronaviruses (SARS-CoV-2)",

    # Rhinoviruses
    "Rhinovirus B37": "Rhinoviruses",
    "Rhinovirus A34": "Rhinoviruses",
    "Rhinovirus C42": "Rhinoviruses",
    "Rhinovirus C2": "Rhinoviruses",
    "Rhinovirus C56": "Rhinoviruses",
    "Rhinovirus A94": "Rhinoviruses",
    "Rhinovirus C36": "Rhinoviruses",

    # Mononegavirales
    "RSVA": "Mononegavirales",
    "RSVB": "Mononegavirales",
    "HPIV4": "Mononegavirales"
}

# Load taxonomy names
taxid_names = {}
with open("index/20250314.taxonomy-names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid not in taxid_names or name_class == "scientific name":
            taxid_names[taxid] = name


# Getting sample pool size
sample_pool_size = {}
with open("[2024] Zephyr sample log - Sampling runs.tsv", "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_date = datetime.strptime(row["date collected"], "%Y-%m-%d")
        if sample_date < datetime(2025, 1, 1) or sample_date > datetime(2025, 2, 25):
            continue
        sample_pool_size[row["sample name"]] = row["total swabs"]


# Parsing virus reads (both those that we didn't validate on purpose and those that we validated)
pooled_data = {}
pathogens = set()

def is_date_in_range(date):
    return datetime(2025, 1, 7) <= date <= datetime(2025, 2, 25)

for sample in sample_pool_size:
    pooled_data[sample] = defaultdict(bool)

for row in csv.DictReader(open(os.path.join(validation_output_dir, "swabs-classified-dedup-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if not is_date_in_range(date):
        continue
    location = row["loc"]
    sample = date.strftime("%y%m%d") + "-" + location + "-NAS"
    pathogen = row["genome_name"]
    pooled_data[sample][pathogen] = True
    pathogens.add(pathogen)

for row in csv.DictReader(open(os.path.join(validation_output_dir, "swabs-non-validated-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if not is_date_in_range(date):
        continue
    location = row["loc"]
    sample = date.strftime("%y%m%d") + "-" + location + "-NAS"
    pathogen = taxid_names[int(row["taxid"])]
    pooled_data[sample][pathogen] = True
    pathogens.add(pathogen)


# Creating virus clade groupings
clades = set(clade_mapping.values())

# Sorting pathogens and clades
pathogens = sorted(list(pathogens))
clades = sorted(list(clades))

# Writing pathogen presence data to file
with open("pathogen_presence.tsv", "w") as f:
    f.write("\t".join(["sample", "pool_size"] + pathogens + clades) + "\n")
    for sample in sorted(pooled_data):
        # Set True if any pathogen in the clade is True
        clade_positivity = defaultdict(bool)
        row = [sample, sample_pool_size[sample]]

        # First determine pathogen status and update clade_positivity
        for pathogen in pathogens:
            status = pooled_data[sample].get(pathogen, False)
            row.append("1" if status else "0")
            clade = clade_mapping[pathogen]
            if status:
                clade_positivity[clade] = True

        # Then add clade statuses to the row
        for clade in clades:
            row.append("1" if clade_positivity[clade] else "0")

        f.write("\t".join(row) + "\n")

