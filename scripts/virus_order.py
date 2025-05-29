#!/usr/bin/env python3
import pandas as pd
import numpy as np

def get_species_order():
    """
    Get species ordering based on group average relative abundance and species median relative abundance.
    
    This function replicates the ordering logic from plot_ra_and_presence.py to ensure
    consistent species ordering across all visualization scripts.
    
    Returns:
        tuple: (panel_df, ordered_species) where:
            - panel_df: DataFrame with species, group, and calculated metrics
            - ordered_species: List of species names in the determined order
    """
    # Load wastewater relative abundance data
    ww = (
        pd.read_csv("tables/ww-ra-summary.tsv", sep="\t")
        .drop(columns="assignment")
        .groupby(["date", "location", "group", "species"], as_index=False)
        .agg({"dedup_hv": "sum", "non_dedup_hv": "sum", "all_reads": "first"})
    )

    # Create complete grid of all date/location/species/group combinations
    grid = ww[["date", "location", "all_reads"]].drop_duplicates().merge(
        ww[["species", "group"]].drop_duplicates(), how="cross"
    )

    # Merge with actual data, filling zeros where no detection occurred
    ww_complete = (
        grid.merge(ww, how="left")
        .fillna({"dedup_hv": 0, "non_dedup_hv": 0})
    )

    # Calculate relative abundance using dedup_hv / all_reads
    ww_complete["relative_abundance"] = ww_complete["dedup_hv"] / ww_complete["all_reads"]

    # Load swab pathogen presence data
    swabs = pd.read_csv("tables/swabs-ra-summary.tsv", sep="\t")
    swabs_df = swabs.drop(columns="assignment").groupby(['date', 'location', 'species', 'group']).agg({
        'dedup_hv': 'sum',
        'non_dedup_hv': 'sum',
        'all_reads': 'first'
    }).reset_index()

    # Create a presence indicator (1 if virus is present)
    swabs_df["present"] = (swabs_df["dedup_hv"] > 0).astype(int)

    # Get all unique species from wastewater data
    all_species = ww_complete[['species', 'group']].drop_duplicates()

    # Count positive swab pools
    swab_counts = (
        swabs_df[swabs_df['present'] == 1]
        .groupby(['species', 'group'])
        .size()
        .reset_index(name='n_positive_pools')
    )

    # Calculate median relative abundance per species
    median_ra = (
        ww_complete.groupby('species')['relative_abundance']
        .median()
        .rename('median_ra')
        .reset_index()
    )

    # Calculate group average relative abundance
    group_avg_ra = (
        ww_complete.groupby('group')['relative_abundance']
        .mean()
        .rename('group_avg_ra')
        .reset_index()
    )

    # Merge all data, starting with all species from wastewater
    panel_df = (all_species
                .merge(swab_counts, how='left', on=['species', 'group'])
                .merge(median_ra, how='left', on='species')
                .merge(group_avg_ra, how='left', on='group')
                .fillna({'n_positive_pools': 0, 'median_ra': 0, 'group_avg_ra': 0})
    )

    # Sort by group average RA (descending) and within each group sort by species median RA (descending)
    panel_df = panel_df.sort_values(['group_avg_ra', 'median_ra'],
                                ascending=[False, False]
                                ).reset_index(drop=True)

    # Create y-position lookup with explicit positions (ensures alignment)
    num_species = len(panel_df)
    spacing = 1.0  # Consistent spacing
    panel_df['y'] = np.arange(num_species) * spacing

    # Extract ordered species list
    ordered_species = panel_df['species'].tolist()

    return panel_df, ordered_species

def get_species_order_filtered(filter_positive=True, max_species=None):
    """
    Get species ordering with optional filtering.
    
    Args:
        filter_positive (bool): If True, only include species with positive swab pools or non-zero median RA
        max_species (int): Maximum number of species to return (top N)
    
    Returns:
        tuple: (panel_df, ordered_species) with optional filtering applied
    """
    panel_df, ordered_species = get_species_order()
    
    if filter_positive:
        # Only include species with non-zero median RA or positive pools
        panel_df = panel_df[(panel_df['median_ra'] > 0) | (panel_df['n_positive_pools'] > 0)]
    
    if max_species is not None and len(panel_df) > max_species:
        panel_df = panel_df.iloc[:max_species]
    
    # Update ordered species list and y positions
    ordered_species = panel_df['species'].tolist()
    panel_df['y'] = np.arange(len(panel_df)) * 1.0
    
    return panel_df, ordered_species

# Color mapping for virus groups (consistent across all scripts)
COLOR_MAPPING = {
    "Coronaviruses (seasonal)": "#05a4a5",
    "Coronaviruses (SARS-CoV-2)": "#445681",
    "Rhinoviruses": "#ba5c97",
    "Mononegavirales": "#8CCEA4",
    "Influenza": "#E08F60",
    "Other": "#9D9D9D"
}

# Desired order for legend display
DESIRED_GROUP_ORDER = [
    "Mononegavirales",
    "Influenza",
    "Coronaviruses (seasonal)",
    "Coronaviruses (SARS-CoV-2)",
    "Rhinoviruses",
]

# Legacy function names for backward compatibility
def get_species_list():
    """Get ordered list of species names."""
    _, ordered_species = get_species_order()
    return ordered_species

def get_species_to_group_mapping():
    """Get mapping from species to group."""
    panel_df, _ = get_species_order()
    return dict(panel_df[['species', 'group']].values)

if __name__ == "__main__":
    # Test the function
    panel_df, ordered_species = get_species_order()
    print(f"Found {len(ordered_species)} species:")
    for i, species in enumerate(ordered_species[:10]):  # Show first 10
        row = panel_df[panel_df['species'] == species].iloc[0]
        print(f"{i+1:2d}. {species:30s} (group: {row['group']:20s}, group_avg_ra: {row['group_avg_ra']:.2e}, median_ra: {row['median_ra']:.2e})")
    
    if len(ordered_species) > 10:
        print(f"... and {len(ordered_species) - 10} more")