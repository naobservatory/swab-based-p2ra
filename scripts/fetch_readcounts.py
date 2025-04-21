#!/usr/bin/env python3

import os
import csv
import gzip
import subprocess
from collections import defaultdict

# Wastewater deliveries
ww_deliveries = [
    "MJ-2025-01-20-a",
    "MJ-2025-01-20-b",
    "MJ-2025-02-12",
    "MJ-2025-03-01",
    "NAO-BCL-2025-03-03",
]

# Swab deliveries
swab_deliveries = [
    "NAO-ONT-20250120-Zephyr8",
    "NAO-ONT-20250127-Zephyr9",
    "NAO-ONT-20250213-Zephyr10",
    "NAO-ONT-20250213-Zephyr10-QC",
    "NAO-ONT-20250220-Zephyr11",
    "NAO-ONT-20250226-Zephyr10-QC2",
    "NAO-ONT-20250313-Zephyr12",
    "NAO-ONT-20250328-Zephyr12b"
]

def fetch_wastewater_readcounts():
    """Fetch and process read counts from wastewater samples."""
    sample_reads = defaultdict(int)

    for delivery in ww_deliveries:
        # Create the directory structure
        os.makedirs(f"deliveries/{delivery}/output/results", exist_ok=True)

        # Use sync with include flag to get only read_counts files
        subprocess.run([
            "aws", "s3", "sync",
            f"s3://nao-mgs-simon/{delivery}/2.8.1/20250314/output/results/",
            f"deliveries/{delivery}/output/results/",
            "--exclude", "*",
            "--include", "read_counts.tsv.gz"
        ])

        # Process the read counts
        try:
            with gzip.open(f"deliveries/{delivery}/output/results/read_counts.tsv.gz", "rt") as inf:
                for row in csv.DictReader(inf, delimiter="\t"):
                    sample_id = row["sample"]
                    sample_reads[sample_id] += int(row["n_read_pairs"])
        except FileNotFoundError:
            print(f"Warning: read_counts file not found for {delivery}")

    # Write the results
    with open("n_reads_per_ww_sample.tsv", "wt") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(["sample", "reads"])
        for sample, reads in sorted(sample_reads.items()):
            writer.writerow([sample, reads])

    return sample_reads

def fetch_swab_readcounts():
    """Fetch and process read counts from swab samples."""
    sample_reads = defaultdict(int)

    for delivery in swab_deliveries:
        # Create the directory structure
        os.makedirs(f"deliveries/{delivery}/output/results", exist_ok=True)

        # Use sync with include flag to get only read_counts files
        subprocess.run([
            "aws", "s3", "sync",
            f"s3://nao-mgs-simon/v2.8.3.2/{delivery}/output/results/",
            f"deliveries/{delivery}/output/results/",
            "--exclude", "*",
            "--include", "read_counts.tsv.gz"
        ])

        # Process the read counts
        try:
            with gzip.open(f"deliveries/{delivery}/output/results/read_counts.tsv.gz", "rt") as inf:
                for row in csv.DictReader(inf, delimiter="\t"):
                    sample_id = row["sample"]
                    sample_reads[sample_id] += int(row["n_reads_single"])
        except FileNotFoundError:
            print(f"Warning: read_counts file not found for {delivery}")

    # Write the results
    with open("n_reads_per_swab_sample.tsv", "wt") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(["sample", "reads"])
        for sample, reads in sorted(sample_reads.items()):
            writer.writerow([sample, reads])

    return sample_reads

def main():
    ww_reads = fetch_wastewater_readcounts()
    swab_reads = fetch_swab_readcounts()

    # Write combined results if needed
    with open("n_reads_per_all_samples.tsv", "wt") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(["sample", "reads", "type"])

        for sample, reads in sorted(ww_reads.items()):
            writer.writerow([sample, reads, "wastewater"])

        for sample, reads in sorted(swab_reads.items()):
            writer.writerow([sample, reads, "swab"])

if __name__ == "__main__":
    main()