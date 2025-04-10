#!/usr/bin/env python3
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json
from datetime import datetime
from dateutil import parser
from collections import defaultdict, Counter

validation_output_dir = "validation-output"

taxid_names = {}
with open("index/20250314.taxonomy-names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid not in taxid_names or name_class == "scientific name":
            taxid_names[taxid] = name

dashboard_dir = os.path.expanduser("~/code/mgs-restricted/dashboard")

with open(os.path.join(dashboard_dir, "metadata_samples.json")) as f:
    metadata_samples = json.load(f)

# ============================================================================
#                                Load prevalence
# ============================================================================
virus_prevalence = {}
with open("virus_prevalence_estimates.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        # Clean up pathogen name
        pathogen = row["Virus"].replace(".", " ")

        # Make viruses match the names in the wastewater data.
        if pathogen == "RSVA":
            pathogen = "RSV-A"
        elif pathogen == "RSVB":
            pathogen = "RSV-B"
        if pathogen == "HPIV4":
            pathogen = "HPIV4b"

        # Special handling for coronavirus groups
        if pathogen == "Coronaviruses  SARS CoV 2 ":
            pathogen = "Coronaviruses (SARS-CoV-2)"
        elif pathogen == "Coronaviruses  seasonal ":
            pathogen = "Coronaviruses (seasonal)"

        # Store prevalence data
        virus_prevalence[pathogen] = {
            "prevalence": float(row["Prevalence"]),
            "lower_ci": float(row["Lower_CI"]),
            "upper_ci": float(row["Upper_CI"])
        }



# ============================================================================
#                               Load WW HV data
# ============================================================================

pathogens = set()

sample_pathogens = defaultdict(Counter)

for row in csv.DictReader(open(os.path.join(validation_output_dir, "ww-classified-dedup-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if date < datetime(2025, 1, 7) or date > datetime(2025, 2, 25):
        continue
    sample = row["sample"]
    pathogen = row["genome_name"]
    pathogens.add(pathogen)
    total_reads = metadata_samples[sample]["reads"]
    sample_pathogens[sample][pathogen] += 1


for row in csv.DictReader(open(os.path.join(validation_output_dir, "ww-non-validated-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if date < datetime(2025, 1, 7) or date > datetime(2025, 2, 25):
        continue

    sample = row["sample_id"]
    taxid = row["taxid"]
    pathogen = taxid_names[int(taxid)]
    pathogens.add(pathogen)
    total_reads = metadata_samples[sample]["reads"]
    sample_pathogens[sample][pathogen] += 1

pathogens = sorted(list(pathogens))
print(sample_pathogens)


# We don't want to include non-respiratory pathogens in our analysis.
pathogens_to_ignore = [
    "Apodemus agrarius picornavirus strain Longwan-Rn37 polyprotein",
    "Coxsackievirus A1",
    "Coxsackievirus A19",
    "Coxsackievirus A22",
    "Coxsackievirus A24",
    "Coxsackievirus A4",
    "Coxsackievirus A6",
    "Coxsackievirus A9",
    "Echovirus E11",
    "Enterovirus A",
    "Enterovirus A71",
    "Enterovirus A76",
    "Enterovirus B",
    "Enterovirus C113",
    "Enterovirus C116",
    "Enterovirus C99",
    "Enterovirus D68",
    "Human mastadenovirus A",
    "Human mastadenovirus B114",
    "Poliovirus 2",
    "Human adenovirus 5",
    "Human enterovirus isolate EV/Human/CMRHP58/CMR/2014",
]

respiratory_pathogens = [p for p in pathogens if p not in pathogens_to_ignore]


# Create virus clade groupings
clade_mapping = {
    # Coronaviruses (seasonal)
    "Human coronavirus OC43": "Coronaviruses (seasonal)",
    "Human coronavirus 229E": "Coronaviruses (seasonal)",
    "Human coronavirus HKU1": "Coronaviruses (seasonal)",
    "Human coronavirus NL63": "Coronaviruses (seasonal)",

    # Coronaviruses (SARS-CoV-2)
    "Severe acute respiratory syndrome coronavirus 2": "Coronaviruses (SARS-CoV-2)",

    # Influenza
    "Flu B": "Influenza",
    "H1N1": "Influenza",
    "H3N2": "Influenza",

    # Rhinoviruses
    "Rhinovirus A": "Rhinoviruses",
    "Rhinovirus A1": "Rhinoviruses",
    "Rhinovirus A12": "Rhinoviruses",
    "Rhinovirus A24": "Rhinoviruses",
    "Rhinovirus A40": "Rhinoviruses",
    "Rhinovirus A54": "Rhinoviruses",
    "Rhinovirus A80": "Rhinoviruses",
    "Rhinovirus B3": "Rhinoviruses",
    "Rhinovirus B103": "Rhinoviruses",
    "Rhinovirus C": "Rhinoviruses",
    "Rhinovirus C1": "Rhinoviruses",
    "Rhinovirus C2": "Rhinoviruses",
    "Rhinovirus C3": "Rhinoviruses",
    "Rhinovirus C4": "Rhinoviruses",
    "Rhinovirus C7": "Rhinoviruses",
    "Rhinovirus C8": "Rhinoviruses",
    "Rhinovirus C11": "Rhinoviruses",
    "Rhinovirus C19": "Rhinoviruses",
    "Rhinovirus C20": "Rhinoviruses",
    "Rhinovirus C34": "Rhinoviruses",
    "Rhinovirus C36": "Rhinoviruses",
    "Rhinovirus C44": "Rhinoviruses",
    "Rhinovirus C55": "Rhinoviruses",
    "Rhinovirus C56": "Rhinoviruses",

    # Mononegavirales
    "RSV-A": "Mononegavirales",
    "RSV-B": "Mononegavirales",
    "HMPV-1": "Mononegavirales",
    "HPIV1": "Mononegavirales",
    "HPIV2": "Mononegavirales",
    "HPIV4b": "Mononegavirales"
}

# ============================================================================
#                               Calculate P2RA
# ============================================================================


with open("ww_p2ra.tsv", "w") as f:
    f.write("\t".join(["sample", "total_reads", "pathogen", "n_reads", "relative_abundance", "median_prevalence", "lower_ci", "upper_ci", "median_incidence", "lower_ci_incidence", "upper_ci_incidence", "ra_1_incidence_median", "ra_1_incidence_lower_ci", "ra_1_incidence_upper_ci"]) + "\n")
    for sample in sample_pathogens:
        clade_counts = defaultdict(int)
        total_reads = metadata_samples[sample]["reads"]

        # Process individual pathogens
        for pathogen in respiratory_pathogens:
            row = [sample, total_reads]
            n_reads = sample_pathogens[sample].get(pathogen, 0)

            # Update clade counts
            clade = clade_mapping[pathogen]
            clade_counts[clade] += n_reads

            relative_abundance = n_reads / total_reads
            prevalence_metrics = virus_prevalence.get(pathogen, None)

            row.append(pathogen)
            row.append(n_reads)
            row.append(relative_abundance)

            if prevalence_metrics:
                median_prevalence = prevalence_metrics["prevalence"]
                lower_ci = prevalence_metrics["lower_ci"]
                upper_ci = prevalence_metrics["upper_ci"]

                shedding_duration = 7/7 if "Rhinovirus" in pathogen else \
                                   6/7 if "coronavirus" in pathogen.lower() else \
                                   5/7  # Mononegavirales

                # Calculate weekly incidence
                median_incidence = median_prevalence / (shedding_duration * (1 - median_prevalence))
                lower_ci_incidence = lower_ci / (shedding_duration * (1 - lower_ci))
                upper_ci_incidence = upper_ci / (shedding_duration * (1 - upper_ci))

                ra_1_incidence_median = relative_abundance * (0.01 / median_incidence)
                ra_1_incidence_lower_ci = relative_abundance * (0.01 / lower_ci_incidence)
                ra_1_incidence_upper_ci = relative_abundance * (0.01 / upper_ci_incidence)

                row.extend([median_prevalence, lower_ci, upper_ci,
                           median_incidence, lower_ci_incidence, upper_ci_incidence,
                           ra_1_incidence_median, ra_1_incidence_lower_ci, ra_1_incidence_upper_ci])
                f.write("\t".join([str(item) for item in row]) + "\n")
            else:
                row.extend(["nan", "nan", "nan", "nan", "nan", "nan", "nan", "nan", "nan"])
                f.write("\t".join([str(item) for item in row]) + "\n")
                print(f"No prevalence metrics for {pathogen}")

        # Process clades
        for clade in clade_counts:
            row = [sample, total_reads]
            row.append(clade)
            row.append(clade_counts[clade])

            relative_abundance = clade_counts[clade] / total_reads
            row.append(relative_abundance)

            prevalence_metrics = virus_prevalence.get(clade, None)

            if prevalence_metrics:
                median_prevalence = prevalence_metrics["prevalence"]
                lower_ci = prevalence_metrics["lower_ci"]
                upper_ci = prevalence_metrics["upper_ci"]

                shedding_duration = 7/7 if "Rhinoviruses" in clade else \
                                   6/7 if "Coronaviruses" in clade else \
                                   5/7 if "Mononegavirales" in clade else 5/7

                if shedding_duration == 5/7 and not ("Mononegavirales" in clade):
                    print(f"Unknown clade: {clade}")

                median_incidence = median_prevalence / (shedding_duration * (1 - median_prevalence))
                lower_ci_incidence = lower_ci / (shedding_duration * (1 - lower_ci))
                upper_ci_incidence = upper_ci / (shedding_duration * (1 - upper_ci))

                ra_1_incidence_median = relative_abundance * (0.01 / median_incidence)
                ra_1_incidence_lower_ci = relative_abundance * (0.01 / lower_ci_incidence)
                ra_1_incidence_upper_ci = relative_abundance * (0.01 / upper_ci_incidence)

                row.extend([median_prevalence, lower_ci, upper_ci,
                           median_incidence, lower_ci_incidence, upper_ci_incidence,
                           ra_1_incidence_median, ra_1_incidence_lower_ci, ra_1_incidence_upper_ci])
                f.write("\t".join([str(item) for item in row]) + "\n")
            else:
                print(f"No prevalence metrics for {clade}")
                row.extend(["nan", "nan", "nan", "nan", "nan", "nan", "nan", "nan", "nan"])
                f.write("\t".join([str(item) for item in row]) + "\n")






