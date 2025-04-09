#!/usr/bin/env python3

import os
import csv
import sys
import glob
import gzip
import hashlib
import subprocess
from collections import defaultdict, Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser

VERBOSE = False

target_deliveries = [
    "NAO-BCL-2025-03-03",
    "MJ-2025-01-20-a",
    "MJ-2025-01-20-b",
    "MJ-2025-03-01",
    "MJ-2025-02-12",
]

genome_names = {
    'ON554083.1': 'Human coronavirus OC43',
    'OR833061.1': 'Human coronavirus NL63',
    'OR266946.1': 'Human coronavirus 229E',
    'ON553964.1': 'Human coronavirus HKU1',

    'MT114540.1': 'Canine coronavirus strain', # probably ignore

    'PP756350.1': 'Coxsackievirus A22',
    'PP756349.1': 'Coxsackievirus A22',
    'PQ889404.1': 'Coxsackievirus A24',
    'AB828290.1': 'Coxsackievirus A19',
    'PP711771.1': 'Coxsackievirus A19',

    'KP290111.1': 'Coxsackievirus A9',
    'OP207966.1': 'Coxsackievirus A9',
    'PP585361.1': 'Coxsackievirus A6',
    'PQ057343.1': 'Coxsackievirus A4',
    'JX174176.1': 'Coxsackievirus A1',
    'PP756363.1': 'Coxsackievirus A1',

    'PP621603.1': 'Enterovirus A',
    'PP548259.1': 'Enterovirus A71',
    'MH118028.1': 'Enterovirus A76',

    'MK815002.1': 'Enterovirus B',

    'PQ119791.1': 'Enterovirus C116',
    'KC344833.1': 'Enterovirus C113',
    'PP756365.1': 'Enterovirus C99',

    'PQ766584.1': 'Enterovirus D68',

    'MH933859.1': 'Human enterovirus isolate EV/Human/CMRHP58/CMR/2014',

    'JX275201.2': 'Poliovirus 2',

    'PV178557.1': 'Rhinovirus A',
    'MN749156.1': 'Rhinovirus A1',
    'OP342739.1': 'Rhinovirus A12',
    'LC699415.1': 'Rhinovirus A24',
    'LC720413.1': 'Rhinovirus A40',
    'MZ540950.1': 'Rhinovirus A54',
    'PP194016.1': 'Rhinovirus A62',
    'OK181467.1': 'Rhinovirus A80',

    'PV178382.1': 'Rhinovirus B3',

    'PP991504.1': 'Rhinovirus C',
    'PV178212.1': 'Rhinovirus C',
    'PV178246.1': 'Rhinovirus C',
    'OM001407.1': 'Rhinovirus C1',
    'PV206816.1': 'Rhinovirus C2',
    'PV178562.1': 'Rhinovirus C3',
    'EF582385.1': 'Rhinovirus C4',
    'PV178659.1': 'Rhinovirus C7',
    'PV178426.1': 'Rhinovirus C8',
    'MZ629171.1': 'Rhinovirus C11',
    'MN369038.1': 'Rhinovirus C19',
    'MZ268714.1': 'Rhinovirus C20',
    'PP194074.1': 'Rhinovirus C28',
    'ON729335.1': 'Rhinovirus C34',
    'OP342692.1': 'Rhinovirus C36',
    'MZ153259.1': 'Rhinovirus C44',
    'OM001443.1': 'Rhinovirus C55',
    'MG950179.1': 'Rhinovirus C56',

    'PV178449.1': 'Rhinovirus B103',

    'OQ969163.1': 'Echovirus E11',

    'PV178204.1': 'RSV-A',
    'PQ899930.1': 'RSV-B',

    'PV260519.1': 'HMPV-1',


    'PV206814.1': 'HPIV1',
    'OQ990770.1': 'HPIV2',
    'PP334259.1': 'HPIV3',
    'KY645962.1': 'HPIV4b',

    'PV333623.1': 'H3N2', # A/Massachusetts/ISC-1252/2025(H3N2) S1
    'PV334560.1': 'H3N2', # A/Texas/ISC-1104/2025(H3N2) S2
    'PV333545.1': 'H3N2', # A/Massachusetts/ISC-1244/2025(H3N2) S3
    'PV271515.1': 'H3N2', # A/Washington/GKISBBBE05023/2025(H3N2) S4
    'PV148462.1': 'H3N2', # A/Minnesota/ISC-1276/2024 S5
    'PV261577.1': 'H3N2', # A/Washington/USAFSAM-15402/202 S6
    'PV249902.1': 'H3N2', # Influenza A virus (A/Washington/WA-UW-93632/2024(H3N2))
    'PV149945.1': 'H3N2', # A/Florida/ISC-1308/2024 S8

    'PV152011.1': 'H1N1', # A/Massachusetts/MA-Broad_BWH-17137/2024 S1
    'PV136994.1': 'H1N1', # A/Colorado/ISC-1555/2024 S2
    'PQ599415.1': 'H1N1', # A/Rhode Island/52/2024 S3
    'PQ839821.1': 'H1N1', # A/New York/172/2024 S4
    'PV100537.1': 'H1N1', # A/Michigan/UM-10061963273/2025 S5
    'PV287334.1': 'H1N1', # A/Massachusetts/ISC-1338/2024 S6
    'PQ598835.1': 'H1N1', # A/Minnesota/107/2024(H1N1) S7
    'PQ260637.1': 'H1N1', # A/West Virginia/48/2024 S8

    'PV150683.1': 'Flu B', # S1
    'PV150682.1': 'Flu B', # S2
    'PV150684.1': 'Flu B', # S3
    'PV150685.1': 'Flu B', # S4
    'PV150686.1': 'Flu B', # S5
    'PV150687.1': 'Flu B', # S6
    'PV150688.1': 'Flu B', # S7
    'PV150689.1': 'Flu B', # S8

    'MN901833.1': 'Human mastadenovirus A',
    'PQ189755.1': 'Human mastadenovirus B114',
    'OR876398.1': 'Human adenovirus 5',

    # Some picornavirus that should be dropped
    'NC_076019.1': 'Apodemus agrarius picornavirus strain Longwan-Rn37 polyprotein',
}

genomes = sorted(genome_names)

with open("to_validate_ww.fasta", "w") as outfasta:
    with open("to_validate_ww.tsv", "w") as outtsv:
        header_written = False
        for delivery in target_deliveries:
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


read_tsv = "to_validate_ww.tsv"
reads_fasta = "to_validate_ww.fasta"

read_dates = {}
read_locs = {}
read_seqs = {}
read_taxids = {}
read_samples = {}

with open(read_tsv, "r") as inf, \
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

db_name = short_hexhash("_".join(genomes))
blastdb = "blastdbs/blast_" + db_name + "_db"
if not any(glob.glob("%s*" % blastdb)):
    all_genomes = []
    # Check if blastdbs directory exists and create it if not
    os.makedirs("blastdbs", exist_ok=True)
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

blast_results = "%s.blast" % reads_fasta

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
print(len(good_primary_hits))

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
            dedup_observations[key] = [] # Note, might we want to have the key be the underlying sample? That way we differentiate between MJ and NAO samples.

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

with open("ww-classified-dedup-reads.tsv", "w") as outf:
    # Write header
    header = ["read_id", "date", "loc", "sample", "genome", "genome_name", "taxid", "start", "end", "bit_score"]
    outf.write("\t".join(header) + "\n")

    # Write data rows
    for row in dedup_entries:
        outf.write("\t".join(str(x) for x in row) + "\n")

with open("ww-classified-reads.fasta", "w") as outf:
    for read_id in sorted(read_dates):
        outf.write(">%s %s %s %s %s %s\n%s\n" % (
            read_id,
            read_dates[read_id],
            read_locs[read_id],
            *good_primary_hits.get(read_id, ("unassigned",
                                        "no-bit-score")),
            read_samples[read_id],
            read_seqs[read_id],
            )
        )

missing = set()
for read_id in read_dates:
    if read_id not in good_primary_hits:
        missing.add(read_id)


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

# Check which Bowtie2 taxids are missing
missing_taxids = Counter()
for read_id in missing:
    missing_taxids[read_taxids[read_id]] += 1

if True:
    for count, taxid in sorted((c,t) for (t,c) in missing_taxids.items()):
        print(count, taxid, taxid_names[taxid])
        # with open("%s.not-accounted-for.fasta" % (taxid), "w") as outf:
        #     for read_id in sorted(missing):
        #         if read_taxids[read_id] != taxid:
        #             continue
        #         outf.write(">%s\n%s\n" % (read_id, read_seqs[read_id]))

if missing and False:
    print("Failed to find genomes for:")
    for read_id in sorted(missing):
        print("%s %s" % (read_dates[read_id], read_locs[read_id]))
        print(read_seqs[read_id])

with open("ww-not-accounted-for.fasta", "w") as outf:
    for read_id in sorted(missing):
        outf.write(">%s\n%s\n" % (read_id, read_seqs[read_id]))

n_identified = len(good_primary_hits)
n_total = len(read_dates)
print()

print("Identified %s of %s (%.0f%%)" % (
    n_identified, n_total, 100*(n_identified /n_total)))
