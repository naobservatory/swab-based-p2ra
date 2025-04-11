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
    "MJ-2025-01-20-a",
    "MJ-2025-01-20-b",
    "MJ-2025-02-12",
    "MJ-2025-03-01",
    "NAO-BCL-2025-03-03",
]

SARS_COV_2_TAXID = 2697049


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

# Even before BLASTing, some viruses we already want to exclude (due to these being GI viruses)
taxids_to_exclude = [
    130309,  # Human mastadenovirus F
    318843,  # Almpiwar virus, spurious hit
    11149,  # Transmissible gastroenteritis virus
    42783,  # Coxsackievirus A22, mostly GI
    1134340,  # Enterovirus C116, mostly GI
    627439,  # Human enteric coronavirus strain 4408
]

taxids_to_include = [
    11157,  # Mononegavirales (RSV, Parainfluenza, HMPV)
    11118,  # Coronaviridae (SARS-CoV-2, OC43, NL63, 229E, HKU1)
    12059,  # Enterovirus (Rhinoviruses)
    11308,  # Orthomyxoviridae (Influenza)
    10508,  # Adenoviridae (but we don't include mastadenovirus F)
]

def descends_from_target(taxid, cache={}):
    if taxid not in cache:
        if taxid in taxids_to_exclude:
            cache[taxid] = False
        elif taxid in taxids_to_include:
            cache[taxid] = True
        elif taxid in [0, 1]:
            # Got up to the root without finding any targets
            cache[taxid] = False
        else:
            cache[taxid] = descends_from_target(parents[taxid])

    return cache[taxid]


retain_taxids = set()
with gzip.open("index/20250314.total-virus-db-annotated.tsv.gz", "rt") as inf:
    for row in csv.DictReader(inf, delimiter="\t"):
        if row["infection_status_human"] != "0":
            retain_taxids.add(int(row["taxid"]))

# We don't add SARS-CoV-2 to the validation set because BLAST is disproprtionately slow on SARS-CoV-2 reads (as there are so many reference genomes for SARS-CoV-2)

for delivery in target_deliveries:
    print(delivery)
    taxid_counts = Counter()
    to_validate = []
    sars_reads = []
    subprocess.run(
        [
            "aws",
            "s3",
            "sync",
            f"s3://nao-mgs-simon/{delivery}/2.8.1/20250314/output/results",
            f"deliveries/{delivery}/output/results",
        ]
    )
    os.makedirs(f"delivery_analyses/{delivery}", exist_ok=True)
    with gzip.open(
        f"deliveries/{delivery}/output/results/virus_hits_filtered.tsv.gz", "rt"
    ) as inf:
        for row in csv.DictReader(inf, delimiter="\t"):
            taxid = int(row["bowtie2_taxid_best"])
            if taxid not in retain_taxids or not descends_from_target(taxid):
                continue
            read_id = row["seq_id"]
            sample_id = row["sample"].rsplit("_", 1)[0]

            sample_metadata = metadata_samples[sample_id]
            date = sample_metadata["date"]
            fine_location = sample_metadata["fine_location"]

            if fine_location not in ("DNI", "DSI"):
                continue

            if parser.parse(date) < datetime(2025, 1, 1):
                continue

            query_seq_fwd = row["query_seq_fwd"]
            query_seq_rev = row["query_seq_rev"]
            query_qual_fwd = row["query_qual_fwd"]
            query_qual_rev = row["query_qual_rev"]

            for seq, qual in zip(
                (query_seq_fwd, query_seq_rev),
                (query_qual_fwd, query_qual_rev),
            ):
                if taxid == SARS_COV_2_TAXID:
                    sars_reads.append(
                        (taxid, seq, qual, date, fine_location, read_id, sample_id)
                    )
                else:
                    to_validate.append(
                        (taxid, seq, qual, date, fine_location, read_id, sample_id)
                    )

            taxid_counts[taxid] += 1

    with gzip.open(
        f"delivery_analyses/{delivery}/to_validate.tsv.gz", "wt"
    ) as outf:
        outf.write(
            "\t".join(
                ("taxid", "sequence", "quality", "date", "loc", "read_id", "sample_id")
            )
            + "\n"
        )
        for record in sorted(to_validate):
            outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % record)

    # To make online BLASTing work, write to_validate sequences to FASTA files, 1000 reads per file
    part_num = 1
    for i in range(0, len(to_validate), 1000):
        chunk = to_validate[i : i + 1000]
        with open(
            f"delivery_analyses/{delivery}/to_validate_part_{part_num}.fasta", "w"
        ) as fasta_out:
            for taxid, seq, qual, date, loc, read_id, sample_id in chunk:
                fasta_out.write(f">{read_id}::{loc}::{date}\n{seq}\n")
        part_num += 1

    # For offline BLASTing, write to_validate sequences to a single FASTA file
    with open(
        f"delivery_analyses/{delivery}/to_validate.fasta",
        "w",
    ) as fasta_overall:
        for taxid, seq, qual, date, loc, read_id, sample_id in to_validate:
            fasta_overall.write(f">{read_id}::{loc}::{date}::{sample_id}\n{seq}\n")

    # Write SARS-CoV-2 reads to a separate FASTA file
    with gzip.open(
        f"delivery_analyses/{delivery}/non_validated.tsv.gz", "wt"
    ) as outf:
        outf.write(
            "\t".join(("taxid", "sequence", "quality", "date", "loc", "read_id", "sample_id"))
            + "\n"
        )
        for record in sorted(sars_reads):
            outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % record)

    for count, taxid in sorted((c, t) for (t, c) in taxid_counts.items()):
        print(count, taxid, taxid_names[taxid])
