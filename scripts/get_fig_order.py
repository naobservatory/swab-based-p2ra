#!/usr/bin/env python3
import pandas as pd
import numpy as np

# Function to extract the species order from the RA script data
def extract_species_order():
    # Load the same data sources as the RA script to get identical ordering
    ww_df = pd.read_csv("tables/ww-ra-summary.tsv", sep="\t")
    ww_df["relative_abundance"] = ww_df["dedup_hv"] / ww_df["all_reads"]
    
    sw_df = pd.read_csv("tables/swabs-ra-summary.tsv", sep="\t")
    sw_df["present"] = (sw_df["dedup_hv"] > 0).astype(int)
    
    # Recreate the same processing steps
    swab_counts = (
        sw_df[sw_df['present'] == 1]
          .groupby(['species', 'group'])
          .size()
          .reset_index(name='n_positive_pools')
    )
    
    all_species = (
        pd.concat([ww_df[['species', 'group']], sw_df[['species', 'group']]])
          .drop_duplicates()
    )
    
    median_ra = (
        ww_df.groupby('species')['relative_abundance']
          .median()
          .rename('median_ra')
          .reset_index()
    )
    
    group_avg_ra = (
        ww_df.groupby('group')['relative_abundance']
          .mean()
          .rename('group_avg_ra')
          .reset_index()
    )
    
    panel_df = (all_species
                .merge(swab_counts, how='left', on=['species', 'group'])
                .merge(median_ra, how='left', on='species')
                .merge(group_avg_ra, how='left', on=['group'])
                .fillna({'n_positive_pools': 0, 'median_ra': 0, 'group_avg_ra': 0})
    )
    
    # Sort in the same way
    panel_df = panel_df.sort_values(['group_avg_ra', 'median_ra'],
                                   ascending=[False, False]
                                  ).reset_index(drop=True)
    
    # Only include species with non-zero median RA or positive pools
    panel_df = panel_df[(panel_df['median_ra'] > 0) | (panel_df['n_positive_pools'] > 0)]
    
    # Limit to top N species
    N = 20  # Same limit as in plot_ra_and_presence.py
    if len(panel_df) > N:
        panel_df = panel_df.iloc[:N]
    
    # Return the species in the exact same order
    ordered_species = panel_df['species'].tolist()
    ordered_groups = panel_df['group'].tolist()
    
    print("Ordered species:")
    for i, (species, group) in enumerate(zip(ordered_species, ordered_groups)):
        print(f"{i+1}. {species} ({group})")
    
    return ordered_species

if __name__ == "__main__":
    extract_species_order()