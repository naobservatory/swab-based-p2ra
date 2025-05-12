#! /usr/bin/env python3
import pandas as pd
import csv
import json
import os
from datetime import datetime
from collections import defaultdict
from metadata_utils import first_level_mapping, second_level_mapping, is_date_in_range

validation_output_dir = "validation-output"

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
        if not is_date_in_range(sample_date):
            continue
        if "SAL" in row["sample name"]:
            continue
        sample_pool_size[row["sample name"]] = row["total swabs"]


# Parsing virus reads (both those that we didn't validate on purpose and those that we validated)
pooled_data = {}
pathogens = set()

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
first_level_pathogens = set(first_level_mapping(pathogen) for pathogen in pathogens)
clades = set(second_level_mapping(pathogen) for pathogen in first_level_pathogens)

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
            first_level_pathogen = first_level_mapping(pathogen)
            clade = second_level_mapping(first_level_pathogen)
            if status:
                clade_positivity[clade] = True

        # Then add clade statuses to the row
        for clade in clades:
            row.append("1" if clade_positivity[clade] else "0")

        f.write("\t".join(row) + "\n")

