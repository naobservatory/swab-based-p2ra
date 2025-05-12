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
    return datetime(2025, 1, 3) <= date <= datetime(2025, 4, 1)

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
        "Human rhinovirus 25": "Rhinovirus A",
        "Human rhinovirus 62": "Rhinovirus A",
        "Rhinovirus A": "Rhinovirus A",
        "Rhinovirus A1": "Rhinovirus A",
        "Rhinovirus A12": "Rhinovirus A",
        "Rhinovirus A24": "Rhinovirus A",
        "Rhinovirus A34": "Rhinovirus A",
        "Rhinovirus A40": "Rhinovirus A",
        "Rhinovirus A54": "Rhinovirus A",
        "Rhinovirus A80": "Rhinovirus A",
        "Rhinovirus A89": "Rhinovirus A",
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
        "Human mastadenovirus F": "Human mastadenovirus F",
        "Human adenovirus 5": "Human adenovirus 5",

        "Human alphaherpesvirus 1": "HSV-1",

        # Papillomaviruses
        "Gammapapillomavirus type 22": "HPV-22",
        "Gammapapillomavirus type 123": "HPV-123",
        "Human papillomavirus type 38": "HPV-38",
        "Human papillomavirus type 50": "HPV-50",
        "Human papillomavirus type 57": "HPV-57",
        "Human papillomavirus type 75": "HPV-75",
        "Human papillomavirus type 109": "HPV-109",
        "Human papillomavirus type 194": "HPV-194",
        "Human papillomavirus type 168": "HPV-168",
        "Human papillomavirus 2": "HPV-2",
        "Human papillomavirus 5": "HPV-5",
        "Human papillomavirus 8": "HPV-8",
        "Human papillomavirus 12": "HPV-12",
        "Human papillomavirus 17": "HPV-17",
        "Human papillomavirus 19": "HPV-19",
        "Human papillomavirus 20": "HPV-20",
        "Human papillomavirus 22": "HPV-22",
        "Human papillomavirus 23": "HPV-23",
        "Human papillomavirus 24": "HPV-24",
        "Human papillomavirus 37": "HPV-37",
        "Human papillomavirus 38": "HPV-38",
        "Human papillomavirus 49": "HPV-49",
        "Human papillomavirus 65": "HPV-65",
        "Human papillomavirus 98": "HPV-98",
        "Human papillomavirus 100": "HPV-100",
        "Human papillomavirus 104": "HPV-104",
        "Human papillomavirus 105": "HPV-105",
        "Human papillomavirus 107": "HPV-107",
        "Human papillomavirus 110": "HPV-110",
        "Human papillomavirus 111": "HPV-111",
        "Human papillomavirus 115": "HPV-115",
        "Human papillomavirus 118": "HPV-118",
        "Human papillomavirus 122": "HPV-122",
        "Human papillomavirus 124": "HPV-124",
        "Human papillomavirus 145": "HPV-145",
        "Human papillomavirus 182": "HPV-182",

        # Polyomaviruses
        "Merkel cell polyomavirus": "Merkel cell polyomavirus",
        "Betapolyomavirus hominis": "Betapolyomavirus hominis",
        "JC polyomavirus": "JC polyomavirus",
        "Human polyomavirus 6": "Human polyomavirus 6",

        # Noroviruses
        "Norovirus GII.4": "Norovirus GII.4",

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
        "Human mastadenovirus F": "Adenoviruses",
        "Human adenovirus 5": "Adenoviruses",

        # Herpesviruses
        "HSV-1": "Herpesviruses",

        # Papillomaviruses
        "HPV-2": "Papillomaviruses",
        "HPV-5": "Papillomaviruses",
        "HPV-8": "Papillomaviruses",
        "HPV-12": "Papillomaviruses",
        "HPV-17": "Papillomaviruses",
        "HPV-19": "Papillomaviruses",
        "HPV-20": "Papillomaviruses",
        "HPV-22": "Papillomaviruses",
        "HPV-23": "Papillomaviruses",
        "HPV-24": "Papillomaviruses",
        "HPV-37": "Papillomaviruses",
        "HPV-38": "Papillomaviruses",
        "HPV-49": "Papillomaviruses",
        "HPV-50": "Papillomaviruses",
        "HPV-57": "Papillomaviruses",
        "HPV-65": "Papillomaviruses",
        "HPV-75": "Papillomaviruses",
        "HPV-98": "Papillomaviruses",
        "HPV-100": "Papillomaviruses",
        "HPV-104": "Papillomaviruses",
        "HPV-105": "Papillomaviruses",
        "HPV-107": "Papillomaviruses",
        "HPV-109": "Papillomaviruses",
        "HPV-110": "Papillomaviruses",
        "HPV-111": "Papillomaviruses",
        "HPV-115": "Papillomaviruses",
        "HPV-118": "Papillomaviruses",
        "HPV-122": "Papillomaviruses",
        "HPV-123": "Papillomaviruses",
        "HPV-124": "Papillomaviruses",
        "HPV-145": "Papillomaviruses",
        "HPV-168": "Papillomaviruses",
        "HPV-182": "Papillomaviruses",
        "HPV-194": "Papillomaviruses",

        # Polyomaviruses
        "Merkel cell polyomavirus": "Polyomaviruses",
        "Betapolyomavirus hominis": "Polyomaviruses",
        "JC polyomavirus": "Polyomaviruses",
        "Human polyomavirus 6": "Polyomaviruses",
    }
    return mapping[assignment]


def pathogens_to_ignore():
    # Dropping viruses that aren't respiratory or skin pathogens.
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