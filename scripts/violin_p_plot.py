#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from fig_utils import (
    get_detected_species_order,
    COLOR_MAPPING,
    SMALL_GROUP_ORDER,
    GROUPS_TO_DROP,
)


from metadata_utils import first_level_mapping, second_level_mapping


def load_and_process_data():
    posteriors = pd.read_csv("tables/posteriors.tsv", sep="\t")

    # Get ordered species from virus_order module
    ordered_species = get_detected_species_order()

    # Filter to only include species present in posteriors
    ordered_species = [
        species
        for species in ordered_species
        if second_level_mapping(species) not in GROUPS_TO_DROP
    ]
    filtered_posteriors = posteriors[posteriors["species"].isin(ordered_species)]

    return filtered_posteriors, ordered_species


# Load data
filtered_posteriors, ordered_species = load_and_process_data()

present_groups = set()
# Create species colors using the color mapping from virus_order
species_colors = {}
for species in ordered_species:
    group = second_level_mapping(species)
    color = COLOR_MAPPING[group]
    present_groups.add(group)
    species_colors[species] = color

# Convert p values to percentages and add log-transformed values
filtered_posteriors["p_percent"] = filtered_posteriors["p"] * 100
filtered_posteriors["log_p"] = np.log10(filtered_posteriors["p_percent"])

# Filter out values outside 2nd-98th percentile range for each species
filtered_posteriors = (
    filtered_posteriors.groupby("species")
    .apply(
        lambda x: x[
            (x["log_p"] >= x["log_p"].quantile(0.02))
            & (x["log_p"] <= x["log_p"].quantile(0.98))
        ]
    )
    .reset_index(drop=True)
)

plt.figure(figsize=(10, 3.9))

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
    bw=0.5,

)

# Add median, Q15 and Q85 lines for each species
for i, species in enumerate(ordered_species):
    species_data = filtered_posteriors[filtered_posteriors["species"] == species]
    median = species_data["log_p"].median()
    q15 = species_data["log_p"].quantile(0.15)
    q85 = species_data["log_p"].quantile(0.85)

    # Plot median line
    ax.plot(
        [median, median],
        [i - 0.4, i + 0.4],
        color="white",
        linewidth=0.8,
        zorder=3,
    )
    # Plot Q15 line
    ax.plot(
        [q15, q15],
        [i - 0.4, i + 0.4],
        color="white",
        linewidth=0.5,
        linestyle="--",
        zorder=3,
    )
    # Plot Q85 line
    ax.plot(
        [q85, q85],
        [i - 0.4, i + 0.4],
        color="white",
        linewidth=0.5,
        linestyle="--",
        zorder=3,
    )

# Add vertical grid lines at log values (for percentages)
for x in [-2, -1, 0, 1, 2]:
    ax.axvline(
        x=x, color="lightgray", linestyle="-", linewidth=0.3, alpha=0.5, zorder=0
    )

# Add horizontal grid lines aligned with each species
for y in range(len(ordered_species)):
    ax.axhline(
        y=y, color="lightgray", linestyle="-", linewidth=0.5, alpha=0.5, zorder=0
    )


# Remove top and right spines
sns.despine(ax=ax, top=True, right=True, left=True)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", length=0)

# Set title and labels
ax.set_xlabel("Prevalence (%)", fontsize=12)
ax.set_xlim(-2.5, 1.3)

ax.set_xticks([-2, -1, 0, 1])
ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%"])
ax.set_ylabel("")

present_groups_ordered = [g for g in SMALL_GROUP_ORDER if g in present_groups]
handles = [
    plt.Line2D([0], [0], color=COLOR_MAPPING[group], lw=4)
    for group in present_groups_ordered
]

fig = plt.gcf()
fig.legend(
    handles,
    present_groups_ordered,
    loc="center",
    bbox_to_anchor=(0.53, 0.08),
    ncol=len(present_groups_ordered) if len(present_groups_ordered) <= 4 else 3,
    frameon=False,
    fontsize=10,
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.25)


# Save the figure
plt.savefig("figures/pathogen_p_violin.png", dpi=300)
plt.savefig("figures/pathogen_p_violin.svg")
