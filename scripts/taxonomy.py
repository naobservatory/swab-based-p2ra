#! /usr/bin/env python3
import csv
import gzip
from collections import defaultdict

def load_taxonomy_names(taxonomy_file="index/20250314.taxonomy-names.dmp"):
    """Load taxonomy names from a taxonomy names dump file."""
    taxid_names = {}
    with open(taxonomy_file) as inf:
        for line in inf:
            taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
            taxid = int(taxid)
            if taxid not in taxid_names or name_class == "scientific name":
                taxid_names[taxid] = name
    return taxid_names


def load_human_infecting_taxids(virus_db_file="index/20250314.total-virus-db-annotated.tsv.gz"):
    """Load taxids of viruses that can infect humans."""
    retain_taxids = set()
    with gzip.open(virus_db_file, "rt") as inf:
        for row in csv.DictReader(inf, delimiter="\t"):
            if row["infection_status_human"] != "0":
                retain_taxids.add(int(row["taxid"]))
    return retain_taxids


def load_taxonomy_tree(taxonomy_file="index/20250314.taxonomy-nodes.dmp"):
    """Load taxonomy tree structure from a taxonomy nodes dump file."""
    parents = {}
    children = defaultdict(set)
    with open(taxonomy_file) as inf:
        for line in inf:
            child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split("\t|\t")
            parents[int(child_taxid)] = int(parent_taxid)
            children[int(parent_taxid)].add(int(child_taxid))
    return parents, children
