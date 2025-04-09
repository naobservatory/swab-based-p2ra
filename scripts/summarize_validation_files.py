#! /usr/bin/env python3
import os
import csv
import gzip

ww_target_deliveries = [
    "NAO-BCL-2025-03-03",
    "MJ-2025-01-20-a",
    "MJ-2025-01-20-b",
    "MJ-2025-03-01",
    "MJ-2025-02-12",
]

swab_target_deliveries = [
    "NAO-ONT-20250120-Zephyr8",
    "NAO-ONT-20250127-Zephyr9",
    "NAO-ONT-20250213-Zephyr10",
    "NAO-ONT-20250213-Zephyr10-QC",
    "NAO-ONT-20250220-Zephyr11",
    "NAO-ONT-20250226-Zephyr10-QC2",
    "NAO-ONT-20250313-Zephyr12",
]

work_dir = "validation-work"
output_dir = "validation-output"
os.makedirs(work_dir, exist_ok=True)


with open(os.path.join(work_dir, "to_validate_ww.fasta"), "w") as outfasta, open(os.path.join(work_dir, "to_validate_ww.tsv"), "w") as outtsv, open(os.path.join(work_dir, "non_validated_ww.tsv"), "w") as out_non_val_tsv:
    header_written = False
    for delivery in ww_target_deliveries:
        try:
            with gzip.open(f"delivery_analyses/{delivery}/to_validate.tsv.gz", "rt") as inf:

                reader = csv.DictReader(inf, delimiter='\t')
                if not header_written:
                    outtsv.write('\t'.join(reader.fieldnames) + '\n')
                    header_written = True
                for row in reader:
                    outtsv.write('\t'.join(row.values()) + '\n')
        except FileNotFoundError:
            print(f"Warning: Could not find to_validate.tsv.gz for {delivery}")

        try:
            with open(f"delivery_analyses/{delivery}/to_validate.fasta", "rt") as inf:
                for line in inf:
                    outfasta.write(line)
        except FileNotFoundError:
            print(f"Warning: Could not find to_validate.fasta for {delivery}")

        try:
            with gzip.open(f"delivery_analyses/{delivery}/non_validated.tsv.gz", "rt") as inf:
                for line in inf:
                    out_non_val_tsv.write(line)
        except FileNotFoundError:
            print(f"Warning: Could not find non_validated.tsv.gz for {delivery}")



with open(os.path.join(work_dir, "to_validate_swabs.fasta"), "w") as outfasta, open(os.path.join(work_dir, "to_validate_swabs.tsv"), "w") as outtsv, open(os.path.join(work_dir, "non_validated_swabs.tsv"), "w") as out_non_val_tsv:
        header_written = False
        for delivery in swab_target_deliveries:
            try:
                with gzip.open(f"delivery_analyses/{delivery}/to_validate.tsv.gz", "rt") as inf:

                    reader = csv.DictReader(inf, delimiter='\t')
                    if not header_written:
                        outtsv.write('\t'.join(reader.fieldnames) + '\n')
                        header_written = True
                    for row in reader:
                        outtsv.write('\t'.join(row.values()) + '\n')
            except FileNotFoundError:
                print(f"Warning: Could not find to_validate.tsv.gz for {delivery}")

            try:
                with open(f"delivery_analyses/{delivery}/to_validate.fasta", "rt") as inf:
                    for line in inf:
                        outfasta.write(line)
            except FileNotFoundError:
                print(f"Warning: Could not find to_validate.fasta for {delivery}")

            try:
                with gzip.open(f"delivery_analyses/{delivery}/non_validated.tsv.gz", "rt") as inf:
                    for line in inf:
                        out_non_val_tsv.write(line)
            except FileNotFoundError:
                print(f"Warning: Could not find non_validated.tsv.gz for {delivery}")