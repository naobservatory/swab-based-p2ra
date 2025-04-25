#! /usr/bin/env python3

import csv
from datetime import datetime
import os
from collections import Counter, defaultdict
import json
from metadata_utils import is_date_in_range, parse_count_table
# Constants and paths
TABLE_DIR = "tables"

with open("deliveries.json") as f:
    deliveries = json.load(f)

target_deliveries = deliveries["swab-deliveries"]


# Load counts
read_counts = parse_count_table("n_reads_per_swab_sample")

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
            if not is_date_in_range(date):
                continue
            date_loc_read_counts[(date, location)] += read_counts.get(sample, 0)


# Getting sample pool size
sample_data = defaultdict(int)
with open("[2024] Zephyr sample log - Sampling runs.tsv", "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_name = row["sample name"]
        sample_date = datetime.strptime(row["date collected"], "%Y-%m-%d")
        if not is_date_in_range(sample_date):
            continue
        sample_pool_size = int(row["total swabs"])
        location = row["sample source"]
        read_number = date_loc_read_counts[(sample_date, location)]
        sample_data[sample_name,sample_date, location, read_number] += sample_pool_size


# Writing metadata
with open(os.path.join(TABLE_DIR, "swab-sample-metadata.tsv"), "wt") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["sample", "date", "location", "pool_size", "all_reads"])
    for sample_info, sample_pool_size in sorted(sample_data.items()):
        sample_name, sample_date, location, read_number = sample_info
        sample_date_str = sample_date.strftime("%y%m%d")
        writer.writerow([sample_name, sample_date_str, location, sample_pool_size, read_number])