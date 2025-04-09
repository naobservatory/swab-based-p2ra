#!/usr/bin/env python3

import json
import glob
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser


N_PRINT = 30

target_deliveries = [
    "NAO-ONT-20250120-Zephyr8",
    "NAO-ONT-20250127-Zephyr9",
    "NAO-ONT-20250213-Zephyr10",
    "NAO-ONT-20250213-Zephyr10-QC",
    "NAO-ONT-20250220-Zephyr11",
    "NAO-ONT-20250226-Zephyr10-QC2",
    "NAO-ONT-20250313-Zephyr12",
]
# First, run scripts/summarize_validation_files.py. After first run of blast-swabs.py, change to validation-work/swabs-not-accounted-for.fasta
limit_to = f"validation-work/to_validate_swabs.fasta"

limit_read_ids = set()
with open(limit_to) as inf:
    for read_id, seq in SimpleFastaParser(inf):
        read_id = read_id.split("::")[0]
        limit_read_ids.add(read_id)

limit_len = len(limit_read_ids)
print(f"Number of reads being analyzed: {len(limit_read_ids)}")

title_bit_score_sums = Counter()

read_seqs = {}
with open(limit_to) as inf:
    for read_id, seq in SimpleFastaParser(inf):
        read_id = read_id.split("::")[0]
        read_seqs[read_id] = seq


for delivery in target_deliveries:
    for blast_file in glob.glob(f"delivery_analyses/{delivery}/*Alignment.json"):
        with open(blast_file) as inf:
            blast_results = json.load(inf)
        for result in blast_results.get("BlastOutput2") or blast_results.get("blastOutput2", []):
            query_title = result["report"]["results"]["search"]["query_title"]
            query_title = query_title.split("::")[0]
            if query_title not in limit_read_ids:
                continue

            for hit in result["report"]["results"]["search"]["hits"]:
                for description in hit["description"]:
                    for hsp in hit["hsps"]:
                        try:
                            accession = description["accession"]

                            title = description["title"]
                            bit_score = hsp["bit_score"]

                            title_bit_score_sums[(accession, title)] += bit_score
                        except KeyError:
                            print(f"No further info about entry with gid {description['id']}: {description}")


# Print top hits as before
print("\n=== TOP HITS ===")
for bit_score_sum, (accession, title) in sorted(
        (b,t) for (t,b) in title_bit_score_sums.items())[-N_PRINT:]:
    print(round(bit_score_sum), accession, title, sep="\t")