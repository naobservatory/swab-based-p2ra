#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from virus_order import (
    get_species_order_filtered,
    get_species_to_group_mapping,
    COLOR_MAPPING,
    GROUP_ORDER,
)


def load_posteriors_data():
    df = pd.read_csv("statistics/2025-05-27-p2ra/posteriors.tsv", sep="\t")
    return df


def load_pathogen_presence():
    df = pd.read_csv("tables/pathogen_presence.tsv", sep="\t")
    pathogen_columns = df.columns[2:]
    positive_species = [col for col in pathogen_columns if df[col].sum() > 0]
    return positive_species


# Load data
posteriors = load_posteriors_data()
positive_species = load_pathogen_presence()

# Filter posteriors to only include species that have positive entries
filtered_species = [
    species
    for species in posteriors["species"].unique()
    if species in positive_species or any(species in ps for ps in positive_species)
]
filtered_posteriors = posteriors[posteriors["species"].isin(filtered_species)]

# Get species order and mapping from virus_order module
ordered_species = get_species_order_filtered()
species_to_group = get_species_to_group_mapping()

# Filter to only include species in our posteriors data
ordered_species = [
    sp for sp in ordered_species if sp in filtered_posteriors["species"].unique()
]

# Create species colors using the color mapping from virus_order
species_colors = {}
for species in ordered_species:
    group = species_to_group.get(species, "Other")
    species_colors[species] = COLOR_MAPPING.get(group, COLOR_MAPPING["Other"])

# Convert p values to percentages and add log-transformed values
filtered_posteriors["p_percent"] = filtered_posteriors["p"] * 100
filtered_posteriors["log_p"] = np.log10(filtered_posteriors["p_percent"])

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
    inner=None,
    linewidth=0.5,
    width=0.8,
    density_norm="width",
    cut=0.1,
)

# Add vertical grid lines at log values (for percentages)
for x in [-3, -2, -1, 0, 1, 2]:
    ax.axvline(
        x=x, color="lightgray", linestyle="-", linewidth=0.3, alpha=0.5, zorder=0
    )

# Add horizontal grid lines aligned with each species
for y in range(len(ordered_species)):
    ax.axhline(
        y=y, color="lightgray", linestyle="-", linewidth=0.5, alpha=0.5, zorder=0
    )

# Create a custom legend for the different virus groups
present_groups = [g for g in GROUP_ORDER if g in species_to_group.values()]
handles = [
    plt.Line2D([0], [0], color=COLOR_MAPPING.get(group, "gray"), lw=4)
    for group in present_groups
]

# Remove top and right spines
sns.despine(ax=ax, top=True, right=True, left=True)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", length=0)

# Set title and labels
ax.set_xlabel("Prevalence (%)", fontsize=12)
ax.set_xticks([-3, -2, -1, 0, 1, 2])
ax.set_xticklabels(["0.001%", "0.01%", "0.1%", "1%", "10%", "100%"])
ax.set_ylabel("")

# Add legend at the bottom
fig = plt.gcf()
fig.legend(
    handles,
    present_groups,
    loc="center",
    bbox_to_anchor=(0.53, 0.03),
    ncol=len(present_groups) if len(present_groups) <= 4 else 3,
    frameon=False,
    fontsize=10,
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Create the output directory if it doesn't exist
os.makedirs("figures", exist_ok=True)

# Save the figure
plt.savefig("figures/pathogen_p_violin.png", dpi=300)
print("Saved to figures/pathogen_p_violin.png")
