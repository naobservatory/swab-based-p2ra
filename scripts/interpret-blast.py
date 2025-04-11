#!/usr/bin/env python3

# Script to show which genomes in online BLAST results has the highest cumulative bit score across all BLASTed reads. This script is used in conjunction with scripts/blast.py, like so:
# 1. Run interpret-blast.py on to_validate_*.fasta to understand which genome in the online BLAST results has the highest cumulative bit score across all reads.
# 2. Add that genome to ww-genomes.json or swabs-genomes.json, depending on which set of reads you're analyzing.
# 3. Run blast.py, returning the reads that did not show a good match to said genome (`*-not-accounted-for.fasta`).
# 4. Run interpret-blast.py on `*-not-accounted-for.fasta`, again identifying the genome that now best matches the remaining reads.
# 5. Again, add genome to `*-genomes.json`, and repeat the loop using `*-not-accounted-for.fasta` until you hit diminishing returns.


import json
import glob
from collections import Counter, defaultdict
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

N_PRINT = 30

parser = argparse.ArgumentParser()
parser.add_argument("--label", type=str, choices=["swabs", "ww"], required=True)
args = parser.parse_args()

# First, run scripts/summarize_validation_files.py on to_validate_{args.label}.fasta. After first run of blast.py, change to validation-work/{args.label}-not-accounted-for.fasta
limit_to = f"validation-work/to_validate_{args.label}.fasta"


if args.label == "swabs":
    target_deliveries = [
        "NAO-ONT-20250120-Zephyr8",
        "NAO-ONT-20250127-Zephyr9",
        "NAO-ONT-20250213-Zephyr10",
        "NAO-ONT-20250213-Zephyr10-QC",
        "NAO-ONT-20250220-Zephyr11",
        "NAO-ONT-20250226-Zephyr10-QC2",
        "NAO-ONT-20250313-Zephyr12",
    ]
elif args.label == "ww":
    target_deliveries = [
        "NAO-BCL-2025-03-03",
        "MJ-2025-01-20-a",
        "MJ-2025-01-20-b",
        "MJ-2025-03-01",
        "MJ-2025-02-12",
    ]


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