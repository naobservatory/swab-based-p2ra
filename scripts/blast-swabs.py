#!/usr/bin/env python3

import os
import csv
import glob
import hashlib
import subprocess
from collections import defaultdict, Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser

output_dir = "validation-output"
work_dir = "validation-work"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(work_dir, exist_ok=True)

genome_names = {
    'MW587042.1': 'Human coronavirus OC43',
    'OR833061.1': 'Human coronavirus NL63',
    'DQ415901.1': 'Human coronavirus HKU1',
    'PQ243243.1': 'Human coronavirus 229E',

    'PV178553.1': 'Rhinovirus A34',
    'OK649392.1': 'Rhinovirus A94',

    'PV178233.1': 'Rhinovirus B37',

    'PV178199.1': 'Rhinovirus C2',
    'MZ670596.1': 'Rhinovirus C36',
    'PQ602524.1': 'Rhinovirus C42',
    'MW969528.1': 'Rhinovirus C56',

    'PQ376588.1': 'RSVA',
    'PV206811.1': 'RSVB',

    'NC_012485.1': 'Human papillomavirus type 109',
    'KR816181.1': 'Human papillomavirus type 194',

    'PQ800221.1': 'HPIV4',
}

genomes = sorted(genome_names)
read_tsv_gz = os.path.join(work_dir, "to_validate_swabs.tsv")
reads_fasta = os.path.join(work_dir, "to_blast_swabs.fasta")

read_dates = {}
read_locs = {}
read_seqs = {}
read_taxids = {}
read_samples = {}

with open(read_tsv_gz, "rt") as inf, \
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

good_primary_hits = {}

# First, identify primary genomes for each read
with open(blast_results) as inf:
    for line in inf:
        (read_id, genome, pct_identity, alignment_length,
         mismatches, gap_openings, read_start, read_end,
         genome_start, genome_end, e_value,
         bit_score) = line.strip().split("\t")

        bit_score = float(bit_score)

        if bit_score < 100:
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
            dedup_observations[key] = []

        genome_start = min(int(genome_start_raw), int(genome_end_raw))
        genome_end = max(int(genome_start_raw), int(genome_end_raw))
        assert genome_start <= genome_end

        is_duplicate = False

        # Check if the new observation is within WIGGLE_ROOM bp of anyexisting observation.
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
                genome_start,
                genome_end,
                bit_score,
            ])
        past_observations[key].append((genome_start, genome_end))


# ============================================================================
#                                 OUTPUT
# ============================================================================

with open(os.path.join(output_dir, "swabs-classified-dedup-reads.tsv"), "w") as outf:
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

with open(os.path.join(work_dir, "swabs-not-accounted-for.fasta"), "w") as outf:
    for read_id in sorted(missing):
        outf.write(">%s\n%s\n" % (read_id, read_seqs[read_id]))

# ============================================================================
#                                 VALIDATION
# ============================================================================


if False:
    genome_counts = Counter()
    for read_id, (read_genome, _) in good_primary_hits.items():
        genome_counts[read_genome] += 1
    for count, genome in sorted((c,g) for (g, c) in genome_counts.items()):
        print(count, genome, genome_names[genome], sep="\t")

    for genome in genomes:
        if genome not in genome_counts:
            print("redundant", genome)

if False:
    genome_counts = Counter()
    genomes_for_name = defaultdict(set)
    for read_id, (read_genome, _) in good_primary_hits.items():
        genome_counts[genome_names[read_genome]] += 1
        genomes_for_name[genome_names[read_genome]].add(read_genome)

    date_locs = sorted(set((read_dates[read_id], read_locs[read_id])
                           for read_id in read_locs))
    date_loc_strs = ["%s-%s" % dl for dl in date_locs]

    print("virus", "total_w_dups", "total_unique",
          *date_loc_strs, sep="\t")
    for count, genome_name in sorted((c,g) for (g, c) in genome_counts.items()):

        date_loc_obs = []
        for date, loc in date_locs:
            obs_read_ids = set()
            for (obs_date,
                 obs_loc,
                 obs_genome), obs_list in dedup_observations.items():
                if (obs_genome not in genomes_for_name[genome_name] or
                    obs_date != date or
                    obs_loc != loc):
                    continue
                for read_id, _, _ in obs_list:
                    obs_read_ids.add(read_id)
            date_loc_obs.append(len(obs_read_ids))

        print(genome_name,
              count,
              sum(date_loc_obs),
              *date_loc_obs, sep="\t")


# Check which Minimap2 taxids are missing
missing_taxids = Counter()
for read_id in missing:
    missing_taxids[read_taxids[read_id]] += 1

if True:
    for count, taxid in sorted((c,t) for (t,c) in missing_taxids.items()):
        print(count, taxid, taxid_names[taxid])

if missing and False:
    print("Failed to find genomes for:")
    for read_id in sorted(missing):
        print("%s %s" % (read_dates[read_id], read_locs[read_id]))
        print(read_seqs[read_id])

n_identified = len(good_primary_hits)
n_total = len(read_dates)

# ============================================================================
#                                 SUMMARY
# ============================================================================

print()
print("Identified %s of %s (%.0f%%)" % (
    n_identified, n_total, 100*(n_identified /n_total)))