from datetime import datetime
import csv

def parse_count_table(table_name):
    read_counts = {}
    fname = f"tables/{table_name}.tsv"
    try:
        with open(fname) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                read_counts[row["sample"]] = int(row["reads"])
    except FileNotFoundError:
        raise Exception(f"{fname} not found. Run fetch_readcounts.py first.")
    return read_counts

def is_date_in_range(date):
    return datetime(2025, 1, 3) <= date <= datetime(2025, 2, 15)

def first_level_mapping(assignment):
    mapping = {
        # Coronaviruses (seasonal)
        "Human coronavirus OC43": "HCoV-OC43",
        "Human coronavirus 229E": "HCoV-229E",
        "Human coronavirus HKU1": "HCoV-HKU1",
        "Human coronavirus NL63": "HCoV-NL63",

        # Coronaviruses (SARS-CoV-2)
        "Severe acute respiratory syndrome coronavirus 2": "SARS-CoV-2",

        # Enteroviruses
        "Coxsackievirus A22": "Coxsackievirus A",
        "Coxsackievirus A24": "Coxsackievirus A",
        "Coxsackievirus A19": "Coxsackievirus A",
        "Coxsackievirus A9": "Coxsackievirus A",
        "Coxsackievirus A6": "Coxsackievirus A",
        "Coxsackievirus A4": "Coxsackievirus A",
        "Coxsackievirus A1": "Coxsackievirus A",
        "Enterovirus A": "Enterovirus A",
        "Enterovirus A71": "Enterovirus A",
        "Enterovirus A76": "Enterovirus A",
        "Enterovirus B": "Enterovirus B",
        "Enterovirus C116": "Enterovirus C",
        "Enterovirus C113": "Enterovirus C",
        "Enterovirus C99": "Enterovirus C",
        "Enterovirus D68": "Enterovirus D",
        "Echovirus E11": "Echovirus E",
        "Human enterovirus isolate EV/Human/CMRHP58/CMR/2014": "Human enterovirus",
        "Poliovirus 2": "Poliovirus",

        # Rhinoviruses
        "Rhinovirus A": "Rhinovirus A",
        "Rhinovirus A1": "Rhinovirus A",
        "Rhinovirus A12": "Rhinovirus A",
        "Rhinovirus A24": "Rhinovirus A",
        "Rhinovirus A34": "Rhinovirus A",
        "Rhinovirus A40": "Rhinovirus A",
        "Rhinovirus A54": "Rhinovirus A",
        "Rhinovirus A80": "Rhinovirus A",
        "Rhinovirus A94": "Rhinovirus A",
        "Rhinovirus B3": "Rhinovirus B",
        "Rhinovirus B37": "Rhinovirus B",
        "Rhinovirus B101": "Rhinovirus B",
        "Rhinovirus B103": "Rhinovirus B",
        "Rhinovirus C": "Rhinovirus C",
        "Rhinovirus C1": "Rhinovirus C",
        "Rhinovirus C2": "Rhinovirus C",
        "Rhinovirus C3": "Rhinovirus C",
        "Rhinovirus C4": "Rhinovirus C",
        "Rhinovirus C7": "Rhinovirus C",
        "Rhinovirus C8": "Rhinovirus C",
        "Rhinovirus C11": "Rhinovirus C",
        "Rhinovirus C19": "Rhinovirus C",
        "Rhinovirus C20": "Rhinovirus C",
        "Rhinovirus C34": "Rhinovirus C",
        "Rhinovirus C36": "Rhinovirus C",
        "Rhinovirus C42": "Rhinovirus C",
        "Rhinovirus C44": "Rhinovirus C",
        "Rhinovirus C55": "Rhinovirus C",
        "Rhinovirus C56": "Rhinovirus C",

        # Mononegavirales
        "RSV-A": "RSV-A",
        "RSV-B": "RSV-B",
        "RSVA": "RSV-A",
        "RSVB": "RSV-B",
        "HMPV-1": "HMPV-1",
        "HPIV1": "HPIV1",
        "HPIV2": "HPIV2",
        "HPIV4b": "HPIV4",
        "HPIV4": "HPIV4",
        # Influenza
        "H3N2": "H3N2",
        "H1N1": "H1N1",
        "Flu B": "Flu B",

        # Adenoviruses
        "Human mastadenovirus A": "Human mastadenovirus A",
        "Human mastadenovirus B114": "Human mastadenovirus B114",
        "Human adenovirus 5": "Human adenovirus 5",

        # Other
        "Apodemus agrarius picornavirus strain Longwan-Rn37 polyprotein": "Apodemus agrarius picornavirus"
    }
    return mapping[assignment]

def second_level_mapping(assignment):
    mapping = {
        # Coronaviruses (seasonal)
        "HCoV-OC43": "Coronaviruses (seasonal)",
        "HCoV-229E": "Coronaviruses (seasonal)",
        "HCoV-HKU1": "Coronaviruses (seasonal)",
        "HCoV-NL63": "Coronaviruses (seasonal)",

        # Coronaviruses (SARS-CoV-2)
        "SARS-CoV-2": "Coronaviruses (SARS-CoV-2)",

        # Enteroviruses
        "Coxsackievirus A": "Enteroviruses",
        "Enterovirus A": "Enteroviruses",
        "Enterovirus B": "Enteroviruses",
        "Enterovirus C": "Enteroviruses",
        "Enterovirus D": "Enteroviruses",
        "Echovirus E": "Enteroviruses",
        "Human enterovirus": "Enteroviruses",
        "Poliovirus": "Enteroviruses",

        # Rhinoviruses
        "Rhinovirus A": "Rhinoviruses",
        "Rhinovirus B": "Rhinoviruses",
        "Rhinovirus C": "Rhinoviruses",

        # Mononegavirales
        "RSV-A": "Mononegavirales",
        "RSV-B": "Mononegavirales",
        "HMPV-1": "Mononegavirales",
        "HPIV1": "Mononegavirales",
        "HPIV2": "Mononegavirales",
        "HPIV4": "Mononegavirales",

        # Influenza
        "H3N2": "Influenza",
        "H1N1": "Influenza",
        "Flu B": "Influenza",

        # Adenoviruses
        "Human mastadenovirus A": "Adenoviruses",
        "Human mastadenovirus B114": "Adenoviruses",
        "Human adenovirus 5": "Adenoviruses",

        # Other
        "Apodemus agrarius picornavirus": "Other"
    }
    return mapping[assignment]


def pathogens_to_ignore():
    # We're restricting our analysis to RNA respiratory pathogens. Our goal is to make a comparison between wastewater and swabs, and our swab protocol includes a DNase step.
    return [
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
        "Human mastadenovirus F",
        "Human mastadenovirus B114",
        "Poliovirus 2",
        "Human adenovirus 5",
        "Human enterovirus isolate EV/Human/CMRHP58/CMR/2014",
    ]