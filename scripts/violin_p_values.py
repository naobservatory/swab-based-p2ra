#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict

# Load posteriors data
def load_posteriors_data():
    df = pd.read_csv("statistics/2025-05-28-p2ra/posteriors.tsv", sep="\t")
    return df

# Load pathogen presence data to filter species that have any positive entries
def load_pathogen_presence():
    df = pd.read_csv("tables/pathogen_presence.tsv", sep="\t")

    # Get all columns excluding the first two (sample, pool_size)
    pathogen_columns = df.columns[2:]

    # Find species with at least one positive result
    positive_species = []
    for col in pathogen_columns:
        if df[col].sum() > 0:
            positive_species.append(col)

    return positive_species

# Define color palettes for different virus groups (matching plot_ra_and_presence.py)
color_mapping = {
    "Coronaviruses (seasonal)": "#05a4a5",
    "Coronaviruses (SARS-CoV-2)": "#445681",
    "Rhinoviruses": "#ba5c97",
    "Mononegavirales": "#8CCEA4",
    "Influenza": "#E08F60",
    "Other": "#9D9D9D"
}

# Define desired order of virus groups (matching plot_ra_and_presence.py)
desired_order = [
    "Rhinoviruses",
    "Influenza",
    "Mononegavirales",
    "Coronaviruses (SARS-CoV-2)",
    "Coronaviruses (seasonal)",
]

# Load data
posteriors = load_posteriors_data()
positive_species = load_pathogen_presence()

# Filter posteriors to only include species that have positive entries in pathogen_presence.tsv
filtered_species = []
for species in posteriors['species'].unique():
    # Handle different naming formats
    if species in positive_species or any(species in ps for ps in positive_species):
        filtered_species.append(species)

filtered_posteriors = posteriors[posteriors['species'].isin(filtered_species)]

# Use the exact same species order as in plot_ra_and_presence.py
def get_ra_and_presence_species_order():
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
                                   ascending=[True, True]
                                  ).reset_index(drop=True)

    # Only include species with non-zero median RA or positive pools
    panel_df = panel_df[(panel_df['median_ra'] > 0) | (panel_df['n_positive_pools'] > 0)]

    # Limit to top N species
    N = 20  # Same limit as in plot_ra_and_presence.py
    if len(panel_df) > N:
        panel_df = panel_df.iloc[:N]

    # Return the species in the exact same order
    return panel_df['species'].tolist()

# Get species order from RA script
ra_ordered_species = get_ra_and_presence_species_order()

# Filter to only include species in our posteriors data
ordered_species = [sp for sp in ra_ordered_species if sp in filtered_posteriors['species'].unique()]

# Assign colors to each species based on its group
species_colors = {}
species_to_group = dict(filtered_posteriors[['species', 'group']].drop_duplicates().values)

for species in ordered_species:
    group = species_to_group[species]
    species_colors[species] = color_mapping.get(group, color_mapping["Other"])

# Convert p values to percentages and add log-transformed values
filtered_posteriors['p_percent'] = filtered_posteriors['p'] * 100  # Convert to percentage
filtered_posteriors['log_p'] = np.log10(filtered_posteriors['p_percent'])  # Log of percentage

# Create a new column with species name and group for better labels
filtered_posteriors['species_with_group'] = filtered_posteriors['species'] + " (" + filtered_posteriors['group'] + ")"

# Create the figure
plt.figure(figsize=(10, 5.4))

# Create the violin plot with log-transformed percentage values
ax = sns.violinplot(
    x="log_p",
    y="species",
    data=filtered_posteriors,
    order=ordered_species,
    hue="species",
    palette=species_colors,
    legend=False,
    inner=None,  # No internal box or points
    linewidth=0.5,
    width=0.8,
    density_norm="width",
    cut=0.1,  # Don't extend the violin past the data
)

# Add vertical grid lines at log values (for percentages)
for x in [-3,-2, -1, 0, 1, 2]:
    ax.axvline(x=x, color='lightgray', linestyle='-', linewidth=0.3, alpha=0.5, zorder=0)

# Add horizontal grid lines aligned with each species
for y in range(len(ordered_species)):
    ax.axhline(y=y, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.5, zorder=0)

# Create group separators and labels
group_positions = {}
current_groups = []

for i, species in enumerate(ordered_species):
    group = species_to_group[species]
    if group not in current_groups:
        current_groups.append(group)
        group_positions[group] = i

# Create a custom legend for the different virus groups
handles = []
labels = []
present_groups = [g for g in desired_order if g in species_to_group.values()]

# Add legend entries for groups
for group in present_groups:
    color = color_mapping.get(group, 'gray')
    handles.append(plt.Line2D([0], [0], color=color, lw=4))
    labels.append(group)

# Remove top and right spines
sns.despine(ax=ax, top=True, right=True, left=True)
ax.tick_params(axis='y', length=0)  # Remove y-axis tick marks
ax.tick_params(axis='x', length=0)  # Remove x-axis tick marks

# Set title and labels
ax.set_xlabel('Prevalence (%)', fontsize=12)

# Set x-ticks with proper percentage labels
ax.set_xticks([-3,-2, -1, 0, 1, 2])
ax.set_xticklabels(['0.001%', '0.01%', '0.1%', '1%', '10%', '100%'])
ax.set_ylabel('')

# Add legend at the bottom
fig = plt.gcf()
fig.legend(handles, labels,
          loc='center',
          bbox_to_anchor=(0.53, 0.03),  # Position in the bottom space
          ncol=len(present_groups) if len(present_groups) <= 4 else 3,
          frameon=False,
          fontsize=10)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)  # Make room for legend

# Create the output directory if it doesn't exist
os.makedirs('figures', exist_ok=True)

# Save the figure
plt.savefig('figures/pathogen_p_violin.png', dpi=300)
print("Saved to figures/pathogen_p_violin.png")