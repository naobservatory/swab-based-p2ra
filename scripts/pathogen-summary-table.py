#!/usr/bin/env python3

import csv
import os
from collections import defaultdict, Counter
from statistics import median
from scipy.stats import gmean

# Constants and paths
TABLE_DIR = "tables"

# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------

ww_data = []
with open(os.path.join(TABLE_DIR, "ww-ra-summary.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        row["dedup_hv"] = int(row["dedup_hv"])
        row["all_reads"] = int(row["all_reads"])
        ww_data.append(row)

pool_data = Counter()
swab_samples = set()
with open(os.path.join(TABLE_DIR, "pathogen_presence.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        swab_samples.add(row["sample"])
        for column_name, value in row.items():
            if column_name in ("sample", "pool_size"):
                continue
            # all other columns represent pathogens
            pathogen = column_name
            pool_data[pathogen] += int(value)


total_pools = len(swab_samples)

ww_sample_reads = {}
with open(os.path.join(TABLE_DIR, "wastewater-sample-metadata.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ww_sample_reads[row["sample"]] = int(row["all_reads"])


# Define virus groups
respiratory_groups = {
    "Influenza",
    "Mononegavirales",
    "Coronaviruses (seasonal)",
    "Coronaviruses (SARS-CoV-2)",
    "Rhinoviruses",
    "Adenoviruses",
    "Enteroviruses",
}

skin_groups = {"Papillomaviruses", "Polyomaviruses"}

# ---------------------------------------------------------------------
#   Data processing
# ---------------------------------------------------------------------

# Extract all unique virus groups and species
virus_groups = set()
respiratory_viruses = set()
skin_viruses = set()

for row in ww_data:
    species = row["species"]
    group = row["group"]
    virus_groups.add(group)

    if group in respiratory_groups:
        respiratory_viruses.add(species)
    elif group in skin_groups:
        skin_viruses.add(species)

virus_groups = sorted(virus_groups)
respiratory_viruses = sorted(respiratory_viruses)
skin_viruses = sorted(skin_viruses)

species_hv = defaultdict(lambda: defaultdict(lambda: {"dedup_hv": 0, "all_reads": 0}))
group_hv = defaultdict(lambda: defaultdict(lambda: {"dedup_hv": 0, "all_reads": 0}))


for sample, all_reads in ww_sample_reads.items():
    for virus in respiratory_viruses:
        species_hv[virus][sample]["dedup_hv"] = 0
        species_hv[virus][sample]["all_reads"] = all_reads

    for virus in skin_viruses:
        species_hv[virus][sample]["dedup_hv"] = 0
        species_hv[virus][sample]["all_reads"] = all_reads

    for group in virus_groups:
        group_hv[group][sample]["dedup_hv"] = 0
        group_hv[group][sample]["all_reads"] = all_reads


for row in ww_data:
    sample = f"{row['date']}-{row['location']}"
    species = row["species"]
    group = row["group"]

    # Group by species
    species_hv[species][sample]["dedup_hv"] += row["dedup_hv"]

    # Group by virus group
    group_hv[group][sample]["dedup_hv"] += row["dedup_hv"]


def calculate_ww_stats(grouped_data):
    ra_values = []
    for counts in grouped_data.values():
        ra_values.append(counts["dedup_hv"] / counts["all_reads"])

    non_zero_ra = [ra for ra in ra_values if ra > 0]
    geomean_ra = 0
    positive_samples = f"{len(non_zero_ra)}/{len(ra_values)}"
    if non_zero_ra:
        geomean_ra = gmean(non_zero_ra)

    return {
        "median": median(ra_values),
        "geomean": geomean_ra,
        "min": min(ra_values),
        "max": max(ra_values),
        "positive_samples": positive_samples,
    }


def format_number(number):
    """Format a number in scientific notation or as 0"""
    if number != 0:
        return f"{number:.2E}"
    return "0"


results = []

# Add stats for all pathogen types
for pathogen in respiratory_viruses + skin_viruses + virus_groups:
    stats = calculate_ww_stats(species_hv.get(pathogen, group_hv.get(pathogen, {})))
    swab_presence = f"{pool_data.get(pathogen, 0)}/{total_pools}"
    ww_presence = stats["positive_samples"]

    results.append(
        {
            "Pathogen": pathogen,
            "Median WW RA": format_number(stats["median"]),
            "Geomean (excl. 0 RA samples)": format_number(stats["geomean"]),
            "Min WW RA": format_number(stats["min"]),
            "Max WW RA": format_number(stats["max"]),
            "WW Presence": ww_presence,
            "Swab Presence": swab_presence,
        }
    )

# ---------------------------------------------------------------------
# Format results for output
# ---------------------------------------------------------------------

final_results = []

section_mapping = {
    "Respiratory Viruses": respiratory_viruses,
    "Skin Viruses": skin_viruses,
    "Virus Groups": virus_groups,
}

for section_name, pathogen_set in section_mapping.items():
    # Add section header
    final_results.append(
        {
            "Pathogen": section_name,
            "Median WW RA": "",
            "Geomean (excl. 0 RA samples)": "",
            "Min WW RA": "",
            "Max WW RA": "",
            "WW Presence": "",
            "Swab Presence": "",
        }
    )

    # Add filtered results
    final_results.extend([r for r in results if r["Pathogen"] in pathogen_set])

# ---------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------

# Define output columns
columns = [
    "Pathogen",
    "Median WW RA",
    "Geomean (excl. 0 RA samples)",
    "Min WW RA",
    "Max WW RA",
    "WW Presence",
    "Swab Presence",
]

# Write results to file
with open(os.path.join(TABLE_DIR, "pathogen-summary-table.tsv"), "w") as outf:
    writer = csv.DictWriter(outf, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(final_results)
