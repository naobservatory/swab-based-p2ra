#!/usr/bin/env python3

# Script to BLAST reads against a set of genomes, separately returning reads that are accounted for with a good match and those that are not. This script is used in conjunction with scripts/interpret-blast.py, like so:
# 1. Run interpret-blast.py on to_validate_*.fasta to understand which genome in the online BLAST results has the highest cumulative bit score acorss all reads.
# 2. Add that genome to ww-genomes.json or swabs-genomes.json, depending on which set of reads you're analyzing.
# 3. Run blast.py, returning the reads that did not show a good match to said genome (`*-not-accounted-for.fasta`).
# 4. Run interpret-blast.py on `*-not-accounted-for.fasta`, again identifying the genome that now best matches the remaining reads.
# 5. Again, add genome to `*-genomes.json`, and repeat the loop using `*-not-accounted-for.fasta` until you hit diminishing returns.

import os
import csv
import glob
import hashlib
import subprocess
import argparse
import json
from collections import defaultdict, Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser

output_dir = "validation-output"
work_dir = "validation-work"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(work_dir, exist_ok=True)

# After first using a score of 200, reads mapping to a viruses added to `*-genomes.json` still ended up in the `*-not-accounted-for.fasta`. Many of these reads were short, but did actually have a best match to said genome witha score under 200. Lowering the threshold to 100 fixed this.
BLAST_SCORE_THRESHOLD = 100

parser = argparse.ArgumentParser()
parser.add_argument("--label", type=str, choices=["swabs", "ww"], required=True)
args = parser.parse_args()
label = args.label

with open(os.path.join(work_dir, "%s-genomes.json" % args.label), "r") as f:
    genome_names = json.load(f)
genomes = sorted(genome_names)

read_dates = {}
read_locs = {}
read_seqs = {}
read_taxids = {}
read_samples = {}

reads_fasta = os.path.join(work_dir, "to_blast_%s.fasta" % label)

with open(os.path.join(work_dir, "to_validate_%s.tsv" % label), "rt") as inf, \
     open(reads_fasta, "w") as outf:
    for row in csv.DictReader(inf, delimiter='\t'):
        outf.write(">%s\n%s\n" % (row["read_id"], row["sequence"]))
        read_dates[row["read_id"]] = row["date"]
        read_locs[row["read_id"]] = row["loc"]
        read_seqs[row["read_id"]] = row["sequence"]
        read_taxids[row["read_id"]] = int(row["taxid"])
        read_samples[row["read_id"]] = row["sample_id"]


taxid_names = {}
with open("index/20250314.taxonomy-names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid not in taxid_names or name_class == "scientific name":
            taxid_names[taxid] = name


def short_hexhash(data, algorithm='sha256'):
    return hashlib.new(algorithm, data.encode('utf-8')).hexdigest()[:16]

os.makedirs("blastdbs", exist_ok=True)
db_name = short_hexhash("_".join(genomes))
blastdb = "blastdbs/blast_" + db_name + "_db"
if not any(glob.glob("%s*" % blastdb)):
    all_genomes = []
    allgenomes_fname = "blastdbs/allgenomes_" + db_name + ".fasta"
    genome_download_dir = "genomes"
    os.makedirs(genome_download_dir, exist_ok=True)
    with open(allgenomes_fname, "w") as outf:
        for genome in genomes:
            genome_fname = os.path.join(genome_download_dir, "%s.fasta" % genome)
            if not os.path.exists(genome_fname):
                print("Downloading %s" % genome)
                subprocess.check_call(["scripts/download_fasta.sh", genome, genome_download_dir])
            with open(genome_fname) as inf:
                for title, seq in SimpleFastaParser(inf):
                    outf.write(">%s\n%s\n" % (title, seq))

    subprocess.check_call([
        "makeblastdb",
        "-in", allgenomes_fname,
        "-dbtype", "nucl",
        "-out", blastdb])



blast_results = os.path.join("%s.blast" % reads_fasta)

subprocess.check_call([
    "blastn",
    "-query", reads_fasta,
    "-db", blastdb,
    "-out", blast_results,
    "-outfmt", "6"])

# First, identify primary genomes for each read

good_primary_hits = {}

with open(blast_results) as inf:
    for line in inf:
        (read_id, genome, pct_identity, alignment_length,
         mismatches, gap_openings, read_start, read_end,
         genome_start, genome_end, e_value,
         bit_score) = line.strip().split("\t")

        bit_score = float(bit_score)

        if bit_score < BLAST_SCORE_THRESHOLD:
            continue

        if (read_id not in good_primary_hits or
            good_primary_hits[read_id][-1] < bit_score):
            good_primary_hits[read_id] = genome, bit_score


dedup_observations = {} # (date, loc, genome) -> [(read_id, start, end)]
past_observations = {} # (date, loc, genome) -> [(start, end)]

# Then count each read towards at most one genome

WIGGLE_ROOM=3 # count dups even if start/end is off by this many bp

dedup_entries = []
with open(blast_results) as inf:
    for line in inf:
        (read_id, genome, pct_identity, alignment_length,
         mismatches, gap_openings, read_start, read_end,
         genome_start_raw, genome_end_raw, e_value,
         bit_score) = line.strip().split("\t")

        # Ignore reads that had no hit with a good alignment
        if read_id not in good_primary_hits:
            continue

        primary_genome, _ = good_primary_hits[read_id]

        # Ignore matches to non-primary genomes.
        if genome != primary_genome:
            continue

        key = read_dates[read_id], read_locs[read_id], genome
        if key not in past_observations:
            past_observations[key] = []
            dedup_observations[key] = [] # Note, might we want to have the key be the underlying sample? That way we differentiate between MU and NAO samples.

        genome_start = min(int(genome_start_raw), int(genome_end_raw))
        genome_end = max(int(genome_start_raw), int(genome_end_raw))
        assert genome_start <= genome_end

        is_duplicate = False

        # Check if the new observation is within WIGGLE_ROOM bp of any existing observation.
        for existing_start, existing_end in past_observations[key]:
            if (abs(existing_start - genome_start) <= WIGGLE_ROOM and
                abs(existing_end - genome_end) <= WIGGLE_ROOM):
                is_duplicate = True

        if not is_duplicate:
            dedup_observations[key].append((
                read_id, genome_start, genome_end))
            dedup_entries.append([
                read_id,
                read_dates[read_id],
                read_locs[read_id],
                read_samples[read_id],
                genome,
                genome_names[genome],
                read_taxids[read_id],
                genome_start,
                genome_end,
                bit_score,
            ])
        past_observations[key].append((genome_start, genome_end))

# ============================================================================
#                                 OUTPUT
# ============================================================================


with open(os.path.join(output_dir, "%s-classified-dedup-reads.tsv" % label), "w") as outf:
    # Write header
    header = ["read_id", "date", "loc", "sample", "genome", "genome_name", "taxid", "start", "end", "bit_score"]
    outf.write("\t".join(header) + "\n")

    # Write data rows
    for row in dedup_entries:
        outf.write("\t".join(str(x) for x in row) + "\n")

missing = set()
for read_id in read_dates:
    if read_id not in good_primary_hits:
        missing.add(read_id)

with open(os.path.join(work_dir, "%s-not-accounted-for.fasta" % label), "w") as outf:
    for read_id in sorted(missing):
        outf.write(">%s\n%s\n" % (read_id, read_seqs[read_id]))


# ============================================================================
#                                 SUMMARY
# ============================================================================

print("Genomes identified:")
genome_counts = Counter()
for read_id, (read_genome, _) in good_primary_hits.items():
    genome_counts[read_genome] += 1
for count, genome in sorted((c,g) for (g, c) in genome_counts.items()):
    print(count, genome, genome_names[genome], sep="\t")

# Check which taxids are missing
missing_taxids = Counter()
for read_id in missing:
    missing_taxids[read_taxids[read_id]] += 1

print("Taxids missing:")
for count, taxid in sorted((c,t) for (t,c) in missing_taxids.items()):
    print(count, taxid, taxid_names[taxid])

n_identified = len(good_primary_hits)
n_total = len(read_dates)


print()
print("Identified %s of %s (%.0f%%)" % (
    n_identified, n_total, 100*(n_identified /n_total)))