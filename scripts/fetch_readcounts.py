#!/usr/bin/env python3

import os
import csv
import gzip
import subprocess
from collections import Counter

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

def fetch_readcounts():
    """Fetch and process read counts from both wastewater and swab samples."""
    ww_reads = Counter()
    swab_reads = Counter()

    # Process wastewater deliveries
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
                    ww_reads[sample_id] += int(row["n_read_pairs"])
        except FileNotFoundError:
            raise ValueError(f"Read counts file not found for {delivery}")

    # Write the wastewater results
    with open("n_reads_per_ww_sample.tsv", "wt") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(["sample", "reads"])
        for sample, reads in sorted(ww_reads.items()):
            writer.writerow([sample, reads])

    # Process swab deliveries
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
                    swab_reads[sample_id] += int(row["n_reads_single"])
        except FileNotFoundError:
            print(f"Warning: read_counts file not found for {delivery}")

    # Write the swab results
    with open("n_reads_per_swab_sample.tsv", "wt") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(["sample", "reads"])
        for sample, reads in sorted(swab_reads.items()):
            writer.writerow([sample, reads])

    return ww_reads, swab_reads

def main():
    ww_reads, swab_reads = fetch_readcounts()

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