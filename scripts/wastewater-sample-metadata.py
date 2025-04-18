#! /usr/bin/env python3

import csv
from datetime import datetime
import os
from collections import Counter
from dateutil import parser

# Constants and paths
tables_dir = "tables"
os.makedirs(tables_dir, exist_ok=True)

target_deliveries = [
    "MJ-2025-01-20-a",
    "MJ-2025-01-20-b",
    "MJ-2025-02-12",
    "MJ-2025-03-01",
    "NAO-BCL-2025-03-03",
]


# Load counts
read_counts = {}
with open(os.path.join("outputs", "n_reads_per_ww_sample.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        read_counts[row["sample"]] = int(row["reads"])

data = set()

date_loc_read_counts = Counter()

# Reading metadata
for delivery in target_deliveries:
    with open(
        os.path.join("..", "mgs-metadata", "deliveries", delivery, "metadata.tsv")
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            fine_location = row["fine_location"]
            sample_id = row["sample"]

            if "BCL" in sample_id:
                continue

            if fine_location not in ("DNI", "DSI"):
                continue

            if parser.parse(row["date"]) < datetime(2025, 1, 3):
                continue

            date = datetime.strptime(row["date"], "%Y-%m-%d")
            date_loc_read_counts[(date, fine_location)] += read_counts.get(sample, 0)
            sample = f"{date.strftime('%y%m%d')}-{fine_location}"
            data.add((sample, date, fine_location))

# Writing metadata
with open(os.path.join(tables_dir, "wastewater-sample-metadata.tsv"), "wt") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["sample", "date", "location", "all_reads"])
    for sample in sorted(data):
        sample_name, sample_date, location = sample
        read_count = date_loc_read_counts[(sample_date, location)]
        sample_date_str = sample_date.strftime("%y%m%d")
        writer.writerow([sample_name, sample_date_str, location, read_count])