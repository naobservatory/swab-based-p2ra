#! /usr/bin/env python3

import csv
from datetime import datetime
import os
from collections import Counter, defaultdict

# Constants and paths
tables_dir = "tables"
os.makedirs(tables_dir, exist_ok=True)
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


# Load counts
read_counts = {}
with open(os.path.join("outputs", "n_reads_per_swab_sample.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        read_counts[row["sample"]] = int(row["reads"])

date_loc_read_counts = Counter()
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


# Getting sample pool size
sample_data = defaultdict(int)
with open("[2024] Zephyr sample log - Sampling runs.tsv", "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_name = row["sample name"]
        sample_date = datetime.strptime(row["date collected"], "%Y-%m-%d")
        if sample_date < datetime(2025, 1, 1) or sample_date > datetime(2025, 2, 25):
            continue
        sample_pool_size = int(row["total swabs"])
        location = row["sample source"]
        read_number = date_loc_read_counts[(sample_date, location)]
        sample_data[sample_name,sample_date, location, read_number] += sample_pool_size


# Writing metadata
with open(os.path.join(tables_dir, "swab-sample-metadata.tsv"), "wt") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["sample", "date", "location", "pool_size", "all_reads"])
    for sample_info, sample_pool_size in sorted(sample_data.items()):
        sample_name, sample_date, location, read_number = sample_info
        sample_date_str = sample_date.strftime("%y%m%d")
        writer.writerow([sample_name, sample_date_str, location, sample_pool_size, read_number])