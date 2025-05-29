#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict
from get_fig_order import extract_species_order

# Load and process data
def load_and_process_data():
    # Load data sources
    posteriors = pd.read_csv("statistics/2025-05-28-p2ra/posteriors.tsv", sep="\t")

    # Get ordered species from get_fig_order.py
    ordered_species = extract_species_order()

    # Filter ordered_species to only include those present in posteriors
    ordered_species = [species for species in ordered_species if species in posteriors['species'].unique()]
    # Filter posteriors to only include ordered species
    filtered_posteriors = posteriors[posteriors['species'].isin(ordered_species)]

    return filtered_posteriors, ordered_species

# Define color mapping for virus groups
color_mapping = {
    "Coronaviruses (seasonal)": "#05a4a5",
    "Coronaviruses (SARS-CoV-2)": "#445681",
    "Rhinoviruses": "#ba5c97",
    "Mononegavirales": "#8CCEA4",
    "Influenza": "#E08F60",
}

# Load and process data
filtered_posteriors, ordered_species = load_and_process_data()

# Get present groups and their colors
species_to_group = dict(filtered_posteriors[['species', 'group']].drop_duplicates().values)
present_groups = list(set(species_to_group.values()))
species_colors = {species: color_mapping[species_to_group[species]]
                 for species in ordered_species}

# Calculate mu_ww * 0.01 and add log-transformed values
filtered_posteriors['scaled_mu_ww'] = filtered_posteriors['mu_ww'] * 0.01
filtered_posteriors['log_scaled_mu_ww'] = np.log10(filtered_posteriors['scaled_mu_ww'])

# Create the figure
plt.figure(figsize=(10, 5.4))

# Create the violin plot
ax = sns.violinplot(
    x="log_scaled_mu_ww",
    y="species",
    data=filtered_posteriors,
    order=ordered_species,
    hue="species",
    palette=species_colors,
    legend=False,
    inner=None,
    linewidth=0.5,
    width=0.8,
    density_norm="width",
    cut=0.1,
)

# Configure plot appearance
ax.set_xlabel('RA(1%)', fontsize=12)
ax.set_ylabel('')
ax.set_xticks([-11, -10, -9, -8, -7, -6, -5, -4, -3])
ax.set_xticklabels(['10⁻¹¹', '10⁻¹⁰', '10⁻⁹', '10⁻⁸', '10⁻⁷', '10⁻⁶', '10⁻⁵', '10⁻⁴', '10⁻³'])

# Add vertical grid lines
for x in [-11, -10, -9, -8, -7, -6, -5, -4, -3]:
    ax.axvline(x=x, color='lightgray', linestyle='-', linewidth=0.3, alpha=0.5, zorder=0)

# Remove spines and ticks
sns.despine(ax=ax, top=True, right=True, left=True)
ax.tick_params(axis='y', length=0)
ax.tick_params(axis='x', length=0)

# Create and add legend
handles = [plt.Line2D([0], [0], color=color_mapping[group], lw=4) for group in present_groups]
fig = plt.gcf()
fig.legend(handles, present_groups,
          loc='center',
          bbox_to_anchor=(0.53, 0.03),
          ncol=len(present_groups) if len(present_groups) <= 4 else 3,
          frameon=False,
          fontsize=10)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Create output directory and save figure
os.makedirs('figures', exist_ok=True)
plt.savefig('figures/pathogen_mu_ww_violin.png', dpi=300)
print("Saved to figures/pathogen_mu_ww_violin.png")