#!/usr/bin/env python3

import os
import csv
import gzip
import json
import subprocess
from collections import defaultdict, Counter
from datetime import datetime
from dateutil import parser
import sys
from taxonomy import load_taxonomy_names, load_human_infecting_taxids, load_taxonomy_tree

parents, children = load_taxonomy_tree()
taxid_names = load_taxonomy_names()
retain_taxids = load_human_infecting_taxids()

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

HCOV_299E_TAXID = 11137
SARS_COV_2_TAXID = 2697049

total_reads_to_validate = 0

metadata_samples = {}
for delivery in target_deliveries:
    with open(
        os.path.join("..", "mgs-metadata", "deliveries", delivery, "metadata.tsv")
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            if row.get("demultiplexed") == "False":
                continue
            metadata_samples[sample] = row

genome_ids_to_exclude = [
    "MN369532.1",  # Vaccinia virus isolate, spurious hit
    "AY037928.1",  # Human endogenous retrovirus, spurious hit
    "KY766069.1",  # Zika virus isolate
]

WIGGLE_ROOM=3 # count dups even if start/end is off by this many bp

past_observations = {}

def check_dup(date, fine_location, genome_id, minimap2_ref_start, minimap2_ref_end):
    key = (date, fine_location, genome_id)
    if key not in past_observations:
        past_observations[key] = []

    is_duplicate = False

    for existing_start, existing_end in past_observations[key]:
        if (abs(existing_start - minimap2_ref_start) <= WIGGLE_ROOM and
            abs(existing_end - minimap2_ref_end) <= WIGGLE_ROOM):
            is_duplicate = True

    past_observations[key].append((minimap2_ref_start, minimap2_ref_end))
    return is_duplicate

sample_reads = defaultdict(int)

for delivery in target_deliveries:
    taxid_counts = Counter()
    to_validate = []
    # One sample has a ~7000 HCOV-299E reads, which makes BLAST choke.
    hcov_299E_reads = []
    # We don't add SARS-CoV-2 to the validation set because BLAST chokes on SARS-CoV-2 reads.
    sars_cov_2_reads = []

    past_observations = {}

    subprocess.run(
        [
            "aws",
            "s3",
            "sync",
            f"s3://nao-mgs-simon/v2.8.3.2{delivery}/output/results",
            f"deliveries/{delivery}/output/results",
        ]
    )
    os.makedirs(f"delivery_analyses/{delivery}", exist_ok=True)

    with gzip.open(f"deliveries/{delivery}/output/results/read_counts.tsv.gz", "rt") as inf:
        for row in csv.DictReader(inf, delimiter="\t"):
            sample_id = row["sample"]
            sample_reads[sample_id] += int(row["n_reads_single"])

    with gzip.open(f"deliveries/{delivery}/output/results/hv.tsv.gz", "rt") as inf:
        for row in csv.DictReader(inf, delimiter="\t"):
            taxid = int(row["minimap2_taxid_primary"])
            genome_id = row["minimap2_genome_id_primary"]

            if taxid not in retain_taxids or genome_id in genome_ids_to_exclude:
                continue

            read_id = row["query_name"]

            sample_id = row["sample"]
            sample_metadata = metadata_samples[sample_id]
            date = sample_metadata["date"]
            fine_location = sample_metadata["fine_location"]

            query_seq = row["query_seq_clean"]
            query_qual = row["query_qual_clean"]

            minimap2_ref_start = int(row["minimap2_ref_start"])
            minimap2_ref_end = int(row["minimap2_ref_end"])

            if delivery == "NAO-ONT-20250220-Zephyr11" and taxid == HCOV_299E_TAXID:
                is_duplicate = check_dup(date, fine_location, genome_id, minimap2_ref_start, minimap2_ref_end)
                hcov_299E_reads.append(
                    (
                        taxid,
                        query_seq,
                        query_qual,
                        date,
                        fine_location,
                        read_id,
                        sample_id,
                        is_duplicate,
                    )
                )

            elif taxid == SARS_COV_2_TAXID:
                is_duplicate = check_dup(date, fine_location, genome_id, minimap2_ref_start, minimap2_ref_end)
                sars_cov_2_reads.append(
                    (
                        taxid,
                        query_seq,
                        query_qual,
                        date,
                        fine_location,
                        read_id,
                        sample_id,
                        is_duplicate,
                    )
                )
            else:
                to_validate.append(
                    (
                        taxid,
                        query_seq,
                        query_qual,
                        date,
                        fine_location,
                        read_id,
                        sample_id,
                    )
                )
                total_reads_to_validate += 1
            taxid_counts[taxid] += 1

    with gzip.open(f"delivery_analyses/{delivery}/to_validate.tsv.gz", "wt") as outf:
        outf.write(
            "\t".join(
                ("taxid", "sequence", "quality", "date", "loc", "read_id", "sample_id")
            )
            + "\n"
        )
        for record in sorted(to_validate):
            outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % record)
    with open(
        f"delivery_analyses/{delivery}/to_validate.fasta",
        "w",
    ) as fasta_overall:
        for taxid, seq, qual, date, loc, read_id, sample_id in to_validate:
            fasta_overall.write(f">{read_id}::{loc}::{date}::{sample_id}\n{seq}\n")

    if hcov_299E_reads or sars_cov_2_reads:
        with gzip.open(
            f"delivery_analyses/{delivery}/non_validated.tsv.gz", "wt"
        ) as outf:
            outf.write(
                        "\t".join(("taxid", "sequence", "quality", "date", "loc", "read_id", "sample_id", "is_duplicate"))
                        + "\n"
                    )
            for taxid, seq, qual, date, loc, read_id, sample_id, is_duplicate in hcov_299E_reads:
                outf.write(
                    "\t".join((str(taxid), seq, qual, date, loc, read_id, sample_id, str(is_duplicate)))
                    + "\n"
                )

            for taxid, seq, qual, date, loc, read_id, sample_id, is_duplicate in sars_cov_2_reads:
                outf.write(
                    "\t".join((str(taxid), seq, qual, date, loc, read_id, sample_id, str(is_duplicate)))
                    + "\n"
                )


    for count, taxid in sorted((c, t) for (t, c) in taxid_counts.items()):
        print(count, taxid, taxid_names[taxid])

print(f"Total reads to validate: {total_reads_to_validate}")

with open("n_reads_per_swab_sample.tsv", "wt") as outf:
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(["sample", "reads"])
    for sample, reads in sorted(sample_reads.items()):
        writer.writerow([sample, reads])
