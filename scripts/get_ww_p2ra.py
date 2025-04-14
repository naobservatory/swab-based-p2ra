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



# Shedding duration constants (in days/week)
RHINOVIRUS_SHEDDING = 7/7
# One study found shedding for 10 days.
# https://pubmed.ncbi.nlm.nih.gov/23490188/#:~:text=negative%20for%20HRV%20or%20infected,with%20hypogammaglobulinaemia%20despite%20immunoglobulin%20replacement
# […] Subjects with symptoms of respiratory tract infection, and their
# household contacts, were screened for HRV by reverse transcription PCR. They
# were followed by serial, self-collected nasal swab specimens until negative
# for HRV or infected by another HRV type. We followed 62 HRV infections in 54
# subjects. The mean (95% CI) duration of HRV shedding was 11.4 (8.2-14.7) days
# in children, 10.1 (7.4-12.9) days in adults, and 40.9 (26.4-55.4) days in
# patients with hypogammaglobulinaemia (p <0.001). The duration of respiratory
# tract symptoms correlated with the duration of virus shedding (p 0.002). [….]
#
#Another showed shorter shedding periods (~5 days):
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4347877/#:~:text=Fifty%20(60%25)%20of%20the%2084%20patients%20completing%20all%20four%20study%20visits%20were%20included%20in%20the%20analysis%20of%20viral%20shedding.%20Eleven%20(22%25)%20were%20positive%20for%20HRV%20on%20day%200%20only%2C%2025%20(50%25)%20on%20days%200%20and%203%2C%2012%20(24%25)%20on%20days%200%2C%203%2C%20and%207.%20Only%20two%20(4%25)%20were%20positive%20for%20HRV%20at%20all%20four%20visits%20(Table%203%20).
#
# Let's go with 7 days.

CORONAVIRUS_SHEDDING = 6/7
# One textbook gives ~6-7 days of infection duration.
# https://rcastoragev2.blob.core.windows.net/34ff49fb351d90e526e2439cc212f10f/PMC7149827.pdf#?page=12

#Children seroconvert to HCoV-OC43 and -229E in the first 5 years of life,
#but symptomatic reinfections occur.224 Clinical manifestations of HCoV
#infections are typical of common colds, with average incubation period 1
#day longer than for HRV, and duration of 6–7 days. Low-grade fever may be
#present in up to 20% of patients, and cough and sore throat occur
#frequently. More serious infections of the lower respiratory tract have
#been documented.227 In addition, HCoVs have been detected in 8% of
#influenza-like illnesses in frail elderly people in the United States.237

# SARS-CoV-2 shedding duration is ~7 days.
# https://www.nature.com/articles/s41579-022-00822-w
#
# Based on the shedding duration of SARS-CoV-2 (~7 days), and SARS-CoV-2
# probably being more severe, let’s go with 6 days.

MONONEGAVIRALES_SHEDDING = 5/7  # Default for other viruses like RSV, HPIV
# Picked shedding duration of 5 days based on RSV information.
#
# Sanofi: https://pro.campus.sanofi/us/respiratory-syncytial-virus/articles/rsv-viral-shedding
# For example, adults will shed RSV for 3 to 7 days following the
# infection, while infants generally shed for up to 14 days in mild
# infections, although this could be as long as 3 weeks in infants
# with severe infection and several months for those who are
# immunocompromised.6

# CDC: https://www.cdc.gov/rsv/causes/index.html#:~:text=People%20with%20RSV%20are%20usually%20contagious%20for%203%20to%208%20days%20and%20may%20become%20contagious%20a%20day%20or%20two%20before%20they%20start%20showing%20signs%20of%20illness.%20However%2C%20some%20infants%20and%20people%20with%20weakened%20immune%20systems%20can%20continue%20to%20spread%20the%20virus%20for%204%20weeks%20or%20longer%2C%20even%20after%20they%20stop%20showing%20symptoms.
# People with RSV are usually contagious for 3 to 8 days and may
#  become contagious a day or two before they start showing signs
#  of illness. However, some infants and people with weakened immune
#  systems can continue to spread the virus for 4 weeks or longer,
# even after they stop showing symptoms.

clade_shedding_duration = {
    "Rhinoviruses": RHINOVIRUS_SHEDDING,
    "Coronaviruses (other)": CORONAVIRUS_SHEDDING,
    "Coronaviruses (SARS-CoV-2)": CORONAVIRUS_SHEDDING,
    "Mononegavirales": MONONEGAVIRALES_SHEDDING
}

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
            pathogen = "Coronaviruses (other)"

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


def is_date_in_range(date):
    return datetime(2025, 1, 7) <= date <= datetime(2025, 2, 25)


for row in csv.DictReader(open(os.path.join(validation_output_dir, "ww-classified-dedup-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if not is_date_in_range(date):
        continue
    sample = row["sample"]
    pathogen = row["genome_name"]
    pathogens.add(pathogen)
    total_reads = metadata_samples[sample]["reads"]
    sample_pathogens[sample][pathogen] += 1


for row in csv.DictReader(open(os.path.join(validation_output_dir, "ww-non-validated-reads.tsv")), delimiter="\t"):
    date = datetime.strptime(row["date"], "%Y-%m-%d")
    if not is_date_in_range(date):
        continue

    sample = row["sample_id"]
    taxid = row["taxid"]
    pathogen = taxid_names[int(taxid)]
    pathogens.add(pathogen)
    total_reads = metadata_samples[sample]["reads"]
    sample_pathogens[sample][pathogen] += 1

pathogens = sorted(list(pathogens))
print(sample_pathogens)


# We're restricting our analysis to respiratory pathogens
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
    "Human coronavirus OC43": "Coronaviruses (other)",
    "Human coronavirus 229E": "Coronaviruses (other)",
    "Human coronavirus HKU1": "Coronaviruses (other)",
    "Human coronavirus NL63": "Coronaviruses (other)",

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


def write_pathogen_data(f, sample, total_reads, pathogen, n_reads, virus_prevalence, is_clade):
    row = [sample, total_reads, pathogen, n_reads]

    relative_abundance = n_reads / total_reads
    row.append(relative_abundance)

    prevalence_metrics = virus_prevalence.get(pathogen, None)

    if prevalence_metrics:
        median_prevalence = prevalence_metrics["prevalence"]
        lower_ci = prevalence_metrics["lower_ci"]
        upper_ci = prevalence_metrics["upper_ci"]

        # Determine shedding duration based on pathogen type
        if is_clade:
            clade = pathogen
        else:
            clade = clade_mapping[pathogen]
        shedding_duration = clade_shedding_duration.get(clade, None)
        if shedding_duration is None:
            print(f"Unknown clade: {pathogen}")


        # Calculate weekly incidence
        median_incidence = median_prevalence / (shedding_duration * (1 - median_prevalence))
        lower_ci_incidence = lower_ci / (shedding_duration * (1 - lower_ci))
        upper_ci_incidence = upper_ci / (shedding_duration * (1 - upper_ci))

        # Calculate RA to incidence ratios
        ra_1_incidence_median = relative_abundance * (0.01 / median_incidence)
        ra_1_incidence_lower_ci = relative_abundance * (0.01 / lower_ci_incidence)
        ra_1_incidence_upper_ci = relative_abundance * (0.01 / upper_ci_incidence)

        row.extend([
            median_prevalence, lower_ci, upper_ci,
            median_incidence, lower_ci_incidence, upper_ci_incidence,
            ra_1_incidence_median, ra_1_incidence_lower_ci, ra_1_incidence_upper_ci
        ])
    else:
        print(f"No prevalence metrics for {pathogen}")
        row.extend(["nan"] * 9)  # Add 9 NaN values for missing metrics

    f.write("\t".join([str(item) for item in row]) + "\n")



with open("ww_p2ra.tsv", "w") as f:
    # Write header
    header = [
        "sample", "total_reads", "pathogen", "n_reads", "relative_abundance",
        "median_prevalence", "lower_ci", "upper_ci", "median_incidence",
        "lower_ci_incidence", "upper_ci_incidence", "ra_1_incidence_median",
        "ra_1_incidence_lower_ci", "ra_1_incidence_upper_ci"
    ]
    f.write("\t".join(header) + "\n")

    for sample in sample_pathogens:
        clade_counts = defaultdict(int)
        total_reads = metadata_samples[sample]["reads"]

        # Process individual pathogens and collect clade counts
        for pathogen in respiratory_pathogens:
            n_reads = sample_pathogens[sample].get(pathogen, 0)

            # Update clade counts
            clade = clade_mapping[pathogen]
            clade_counts[clade] += n_reads

            # Write pathogen data
            is_clade = False
            write_pathogen_data(f, sample, total_reads, pathogen, n_reads, virus_prevalence, is_clade)

        # Process clades
        for clade in clade_counts:
            is_clade = True
            write_pathogen_data(f, sample, total_reads, clade, clade_counts[clade], virus_prevalence, is_clade)

