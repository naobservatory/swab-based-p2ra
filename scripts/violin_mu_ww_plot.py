#!/usr/bin/env python3
import json
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
from metadata_utils import second_level_mapping

outside_groups = {
    "Influenza A\n(previously)": "Influenza",
    "SARS-CoV-2\n(previously)": "Coronaviruses (SARS-CoV-2)",
}

outside_colors = {
    "Influenza A\n(previously)": "#ba5c97",
    "SARS-CoV-2\n(previously)": "#05a4a5",
}


def load_posteriors():
    posteriors = pd.read_csv("tables/posteriors.tsv", sep="\t")

    # Get ordered species from virus_order module
    species_order = get_detected_species_order()

    filtered_species = [
        species
        for species in species_order
        if second_level_mapping(species) not in GROUPS_TO_DROP
    ]
    filtered_posteriors = posteriors[posteriors["species"].isin(filtered_species)]
    filtered_posteriors["source"] = "swab-p2ra"

    with open("tables/ww-rai1pct.json", "r") as f:
        previous_estimates = json.load(f)

    # Parse previous estimates for Rothman-2697049 and MU-11320
    outside_species = {
        "Rothman-2697049": "SARS-CoV-2\n(previously)",
        "MU-11320": "Influenza A\n(previously)",
    }
    outside_distributions = {}

    for species, pretty_name in outside_species.items():
        if species in previous_estimates:
            distribution = []
            for key, count in previous_estimates[species].items():
                value = float(key.replace("e", "E"))
                distribution.extend([value] * count)
            outside_distributions[pretty_name] = distribution

    for species, distribution in outside_distributions.items():
        new_rows = pd.DataFrame(
            {
                "species": [species] * len(distribution),
                "scaled_mu_ww": distribution,
                "log_scaled_mu_ww": np.log10(distribution),
                "source": ["previous-estimate"] * len(distribution),
            }
        )

        filtered_posteriors = pd.concat(
            [filtered_posteriors, new_rows], ignore_index=True
        )
        filtered_species.append(species)

    return filtered_posteriors, filtered_species


# Load and process data
filtered_posteriors, ordered_species = load_posteriors()
# Get species to group mapping from virus_order module
groups = set()
species_colors = {}
for species in ordered_species:
    group = outside_groups.get(species)
    color = outside_colors.get(species)
    if group is None:
        group = second_level_mapping(species)
        color = COLOR_MAPPING[group]

    groups.add(group)
    species_colors[species] = color

# Scale mu_ww * 0.01 and add log-transformed values
mask = filtered_posteriors["source"] == "swab-p2ra"
filtered_posteriors.loc[mask, "scaled_mu_ww"] = (
    filtered_posteriors.loc[mask, "mu_ww"] * 0.01
)
filtered_posteriors.loc[mask, "log_scaled_mu_ww"] = np.log10(
    filtered_posteriors.loc[mask, "scaled_mu_ww"]
)

# Filter out values outside 2nd-98th percentile range for each species (for violin styling)
filtered_posteriors = (
    filtered_posteriors.groupby("species")
    .apply(
        lambda x: x[
            (x["log_scaled_mu_ww"] >= x["log_scaled_mu_ww"].quantile(0.02))
            & (x["log_scaled_mu_ww"] <= x["log_scaled_mu_ww"].quantile(0.98))
        ]
    )
    .reset_index(drop=True)
)

# Create the figure
plt.figure(figsize=(10, 6.5))

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
    width=0.6,
    bw=0.5,
)
# Add diagonal shading for the last two violins
for i, species in enumerate(ordered_species[-2:]):
    # Get the violin path
    violin = ax.collections[i + len(ordered_species) - 2]
    path = violin.get_paths()[0]
    vertices = path.vertices

    # Get bounds of violin
    x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()

    # Find points where the diagonal line intersects with the violin
    mask = (vertices[:, 0] >= x_min) & (vertices[:, 0] <= x_max)
    if np.any(mask):
        shade_vertices = vertices[mask]
        x_values = shade_vertices[:, 0]
        y_values = shade_vertices[:, 1]

        # Plot shading
        ax.fill_betweenx(
            y_values,
            x_values,
            x_max,
            color="white",
            alpha=0.3,
            zorder=2,
            edgecolor="none",
        )

# Add horizontal dashed line above second lowest y value
ax.axhline(
    y=len(ordered_species) - 2.47,
    color="black",
    linestyle="--",
    linewidth=0.5,
    alpha=0.5,
    zorder=0,
)


# Add median, Q15 and Q85 lines for each species
for i, species in enumerate(ordered_species):
    species_data = filtered_posteriors[filtered_posteriors["species"] == species]
    median = species_data["log_scaled_mu_ww"].median()
    q15 = species_data["log_scaled_mu_ww"].quantile(0.15)
    q85 = species_data["log_scaled_mu_ww"].quantile(0.85)

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

# Configure plot appearance
ax.set_xlabel("RA(1%)", fontsize=12)
ax.set_ylabel("")
ax.set_xticks([-10, -9, -8, -7, -6, -5, -4])
ax.set_xticklabels(["10⁻¹⁰", "10⁻⁹", "10⁻⁸", "10⁻⁷", "10⁻⁶", "10⁻⁵", "10⁻⁴"])

# Add vertical grid lines
for x in [-10, -9, -8, -7, -6, -5, -4]:
    ax.axvline(
        x=x, color="lightgray", linestyle="-", linewidth=0.3, alpha=0.5, zorder=0
    )

# Remove spines and ticks
sns.despine(ax=ax, top=True, right=True, left=True)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", length=0)

# Create and add legend using GROUP_ORDER
present_groups_ordered = [g for g in SMALL_GROUP_ORDER if g in groups] + ["Influenza"]

handles = [
    plt.Line2D([0], [0], color=COLOR_MAPPING[group], lw=4)
    for group in present_groups_ordered
]

fig = plt.gcf()
fig.legend(
    handles,
    present_groups_ordered,
    loc="center",
    bbox_to_anchor=(0.53, 0.12),
    ncol=len(present_groups_ordered) if len(present_groups_ordered) <= 4 else 3,
    frameon=False,
    fontsize=10,
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.25)

os.makedirs("figures", exist_ok=True)
plt.savefig("figures/pathogen_mu_ww_violin.png", dpi=300)
plt.savefig("figures/pathogen_mu_ww_violin.svg")

