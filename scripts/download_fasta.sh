#!/usr/bin/env bash

set -e
set -u

ACCESSION="$1"
LOCATION="$2"
datasets download virus genome accession "$ACCESSION"
unzip -p ncbi_dataset.zip ncbi_dataset/data/genomic.fna> $LOCATION/$ACCESSION.fasta
rm ncbi_dataset.zip