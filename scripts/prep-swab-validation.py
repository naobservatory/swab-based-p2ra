#!/usr/bin/env python3

import os
import csv
import gzip
import json
import subprocess
from collections import defaultdict, Counter
from datetime import datetime
from dateutil import parser

target_deliveries = [
    "NAO-ONT-20250120-Zephyr8",
    "NAO-ONT-20250127-Zephyr9",
    "NAO-ONT-20250213-Zephyr10",
    "NAO-ONT-20250213-Zephyr10-QC",
    "NAO-ONT-20250220-Zephyr11",
    "NAO-ONT-20250226-Zephyr10-QC2",
    "NAO-ONT-20250313-Zephyr12",
]

dashboard_dir = os.path.expanduser("~/code/mgs-restricted/dashboard")

with open(os.path.join(dashboard_dir, "metadata_samples.json")) as f:
    metadata_samples = json.load(f)

parents = {}
children = defaultdict(set)
with open("index/20250314.taxonomy-nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split("\t|\t")
        parents[int(child_taxid)] = int(parent_taxid)
        children[int(parent_taxid)].add(int(child_taxid))

taxid_names = {}
with open("index/20250314.taxonomy-names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid not in taxid_names or name_class == "scientific name":
            taxid_names[taxid] = name


retain_taxids = set()
with gzip.open("index/20250314.total-virus-db-annotated.tsv.gz", "rt") as inf:
    for row in csv.DictReader(inf, delimiter="\t"):
        if row["infection_status_human"] != "0":
            retain_taxids.add(int(row["taxid"]))

genome_ids_to_exclude = [
    "MN369532.1",  # Vaccinia virus isolate, spurious hit
    "AY037928.1",  # Human endogenous retrovirus, spurious hit
    "KY766069.1",  # Zika virus isolate
]

for delivery in target_deliveries:
    print(delivery)
    taxid_counts = Counter()
    to_validate = []
    # One sample has a ~7000 HCOV-299E reads, which makes BLAST choke.
    hcov_299E_reads = []
    # We don't add SARS-CoV-2 to the validation set because BLAST chokes on SARS-CoV-2 reads.
    sars_cov_2_reads = []
    subprocess.run(
        [
            "aws",
            "s3",
            "sync",
            f"s3://nao-mgs-simon/v2.8.1.3-dev/{delivery}/output/results",
            f"deliveries/{delivery}/output/results",
        ]
    )
    os.makedirs(f"delivery_analyses/{delivery}", exist_ok=True)

    with gzip.open(f"deliveries/{delivery}/output/results/hv.tsv.gz", "rt") as inf:
        for row in csv.DictReader(inf, delimiter="\t"):
            taxid = int(row["minimap2_taxid_primary"])
            genome_id = row["minimap2_genome_id_primary"]

            if taxid not in retain_taxids or genome_id in genome_ids_to_exclude:
                continue

            read_id = row["query_name"]

            sample_id = row["sample"].split("-div")[0]
            sample_metadata = metadata_samples[sample_id]
            date = sample_metadata["date"]
            fine_location = sample_metadata["fine_location"]

            if parser.parse(date) < datetime(2025, 1, 1):
                continue

            query_seq = row["query_seq_clean"]
            query_qual = row["query_qual_clean"]

            if delivery == "NAO-ONT-20250220-Zephyr11" and taxid == 11137:
                hcov_299E_reads.append(
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
            elif taxid == 2697049:
                sars_cov_2_reads.append(
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

    if delivery == "NAO-ONT-20250220-Zephyr11":
        with open(
            f"delivery_analyses/{delivery}/hcov_299E.fasta",
            "w",
        ) as fasta_299E:
            for taxid, seq, qual, date, loc, read_id, sample_id in hcov_299E_reads:
                fasta_299E.write(f">{read_id}::{loc}::{date}::{sample_id}\n{seq}\n")
    if sars_cov_2_reads:
        with open(
            f"delivery_analyses/{delivery}/sars_cov_2.fasta",
            "w",
        ) as fasta_sars_cov_2:
            for taxid, seq, qual, date, loc, read_id, sample_id in sars_cov_2_reads:
                fasta_sars_cov_2.write(
                    f">{read_id}::{loc}::{date}::{sample_id}\n{seq}\n"
                )
    for count, taxid in sorted((c, t) for (t, c) in taxid_counts.items()):
        print(count, taxid, taxid_names[taxid])
