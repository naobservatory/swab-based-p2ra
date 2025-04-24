#! /usr/bin/env python3

import csv
import json
from datetime import datetime
import os
from collections import Counter
from dateutil import parser
from metadata_utils import is_date_in_range

# Constants and paths
TABLE_DIR = "tables"

with open("deliveries.json") as f:
    deliveries = json.load(f)

target_deliveries = deliveries["ww-deliveries"]


# Load counts
read_counts = {}
try:
    with open(os.path.join(TABLE_DIR, "n_reads_per_ww_sample.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            read_counts[row["sample"]] = int(row["reads"])
except FileNotFoundError:
    raise Exception("n_reads_per_ww_sample.tsv not found. Run fetch_readcounts.py first.")

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

            if fine_location not in ("DNI", "DSI"):
                continue

            if not is_date_in_range(parser.parse(row["date"])):
                continue

            date = datetime.strptime(row["date"], "%Y-%m-%d")
            date_loc_read_counts[(date, fine_location)] += read_counts.get(sample, 0)
            sample = f"{date.strftime('%y%m%d')}-{fine_location}"
            data.add((sample, date, fine_location))

# Writing metadata
with open(os.path.join(TABLE_DIR, "wastewater-sample-metadata.tsv"), "wt") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["sample", "date", "location", "all_reads"])
    for sample in sorted(data):
        sample_name, sample_date, location = sample
        read_count = date_loc_read_counts[(sample_date, location)]
        sample_date_str = sample_date.strftime("%y%m%d")
        writer.writerow([sample_name, sample_date_str, location, read_count])