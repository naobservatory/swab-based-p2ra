#!/usr/bin/env python3
import pandas as pd


def get_species_order():
    """Get species ordering: groups by median RA, within groups by mean RA."""
    # Load data
    ww = pd.read_csv("tables/ww-ra-summary.tsv", sep="\t")

    # Calculate relative abundance
    ww = (
        ww.drop(columns="assignment")
        .groupby(["date", "location", "group", "species"], as_index=False)
        .agg({"dedup_hv": "sum", "non_dedup_hv": "sum", "all_reads": "first"})
    )
    ww["relative_abundance"] = ww["dedup_hv"] / ww["all_reads"]

    # Get species info with group
    species_info = ww[["species", "group"]].drop_duplicates()

    # Calculate group median RA
    group_median_ra = (
        ww.groupby("group")["relative_abundance"]
        .median()
        .reset_index()
        .rename(columns={"relative_abundance": "group_median_ra"})
    )

    # Calculate species mean RA
    species_mean_ra = (
        ww.groupby("species")["relative_abundance"]
        .mean()
        .reset_index()
        .rename(columns={"relative_abundance": "species_mean_ra"})
    )

    # Merge and sort
    panel_df = species_info.merge(group_median_ra, on="group").merge(
        species_mean_ra, on="species"
    )

    panel_df = panel_df.sort_values(
        ["group_median_ra", "species_mean_ra"], ascending=[True, True]
    ).reset_index(drop=True)

    species_order = panel_df["species"].tolist()
    return species_order


def get_species_order_filtered():
    """Get filtered species order with only positive detections."""
    # Load swab data to check for positive detections
    swabs = pd.read_csv("tables/swabs-ra-summary.tsv", sep="\t")
    swabs = (
        swabs.drop(columns="assignment")
        .groupby(["species"])
        .agg({"dedup_hv": "sum"})
        .reset_index()
    )
    positive_species = swabs[swabs["dedup_hv"] > 0]["species"].tolist()

    # Get full ordering and filter to positive species
    ordered_species = get_species_order()
    filtered_species = [sp for sp in ordered_species if sp in positive_species]
    return filtered_species


def get_species_to_group_mapping():
    """Get mapping from species to group."""
    ww = pd.read_csv("tables/ww-ra-summary.tsv", sep="\t")
    species_groups = ww[["species", "group"]].drop_duplicates()
    return dict(species_groups.values)


# Color mapping for virus groups
COLOR_MAPPING = {
    "Coronaviruses (seasonal)": "#05a4a5",
    "Coronaviruses (SARS-CoV-2)": "#445681",
    "Rhinoviruses": "#ba5c97",
    "Mononegavirales": "#8CCEA4",
    "Influenza": "#E08F60",
    "Other": "#9D9D9D",
}

# Desired order for legend display
GROUP_ORDER = [
    "Mononegavirales",
    "Influenza",
    "Coronaviruses (seasonal)",
    "Coronaviruses (SARS-CoV-2)",
    "Rhinoviruses",
]
