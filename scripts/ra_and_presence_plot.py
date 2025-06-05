#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import textwrap
from metadata_utils import first_level_mapping, second_level_mapping
from fig_utils import (
    get_species_order,
    COLOR_MAPPING,
    GROUP_ORDER,
)

ZERO_VALUE_REPLACEMENT = 3.16e-10


# Load wastewater relative abundance data
def load_ww_abundance_data():
    df = (
        pd.read_csv("tables/ww-ra-summary.tsv", sep="\t")
        .drop(columns="assignment")
        .groupby(["date", "location", "group", "species"], as_index=False)
        .agg({"dedup_hv": "sum", "non_dedup_hv": "sum", "all_reads": "first"})
    )

    # Create complete grid of all date/location/species/group combinations
    grid = (
        df[["date", "location", "all_reads"]]
        .drop_duplicates()
        .merge(df[["species", "group"]].drop_duplicates(), how="cross")
    )

    # Merge with actual data, filling zeros where no detection occurred
    df_complete = grid.merge(df, how="left").fillna({"dedup_hv": 0, "non_dedup_hv": 0})

    # Calculate relative abundance using dedup_hv / all_reads
    df_complete["relative_abundance"] = (
        df_complete["dedup_hv"] / df_complete["all_reads"]
    )

    return df_complete


# Load swab pathogen presence data
def load_swab_presence_data():
    # Load swab data
    swabs = pd.read_csv("tables/swabs-ra-summary.tsv", sep="\t")
    # Drop assignment column and collapse by summing dedup_hv and all_reads
    swabs_df = (
        swabs.drop(columns="assignment")
        .groupby(["date", "location", "species", "group"])
        .agg({"dedup_hv": "sum", "non_dedup_hv": "sum", "all_reads": "first"})
        .reset_index()
    )

    # Create a presence indicator (1 if virus is present)
    swabs_df["present"] = (swabs_df["dedup_hv"] > 0).astype(int)

    return swabs_df


def handle_zero_values(x_values, y_values, x_scale=0.1, y_scale=0.2):
    x_values = x_values.to_numpy(float)
    mask = x_values == 0
    x_values[mask] = ZERO_VALUE_REPLACEMENT

    x_jitter_multiplier = np.random.uniform(1 - x_scale, 1 + x_scale, mask.sum())
    y_jitter_addition = np.random.uniform(-y_scale, y_scale, mask.sum())
    x_values[mask] *= x_jitter_multiplier
    y_values[mask] += y_jitter_addition

    return x_values, y_values


def main():
    # ---------- Load & prepare data ---------- #
    ww_data = load_ww_abundance_data()
    sw_data = load_swab_presence_data()

    # Get species ordering from virus_order module
    species_to_group = {
        species: second_level_mapping(species) for species in get_species_order()
    }
    ordered_species = get_species_order()[::-1]  # Invert the order of species

    # Create panel dataframe from ordered species list
    panel_df = pd.DataFrame({"species": ordered_species})
    panel_df["group"] = panel_df["species"].map(species_to_group)

    # Count positive swab pools
    swab_counts = (
        sw_data[sw_data["present"] == 1]
        .groupby(["species", "group"])
        .size()
        .reset_index(name="n_positive_pools")
    )

    # Merge swab counts
    panel_df = panel_df.merge(swab_counts, how="left", on=["species", "group"])
    panel_df["n_positive_pools"] = panel_df["n_positive_pools"].fillna(0)

    # Create y-position lookup with explicit positions (to align y positions of both plots)
    num_species = len(panel_df)
    spacing = 1.0
    panel_df["y"] = np.arange(num_species) * spacing

    # ------------ Create figure ---------- #
    fig = plt.figure(figsize=(14, 7.2), dpi=450)
    gs = fig.add_gridspec(
        2, 2, height_ratios=[0.85, 0.15], width_ratios=[3, 1], wspace=0.05, hspace=0
    )
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])

    # Define marker cycling for different groups
    markers = ["s", "o", "^"]  # square, triangle, circle

    # Assign markers to groups
    MARKER_MAPPING = {}
    for i, group in enumerate(panel_df["group"].unique()):
        MARKER_MAPPING[group] = markers[i % len(markers)]

    # ------------ Wastewater plot ---------- #
    for _, row in panel_df.iterrows():
        species = row["species"]
        group = row["group"]
        y_pos = row["y"]
        species_data = ww_data[ww_data["species"] == species]

        # Handle zero values and apply jitter
        x_values = species_data["relative_abundance"]
        y_values = np.full(len(species_data), y_pos)
        x_values, y_values = handle_zero_values(x_values, y_values)

        ax_left.scatter(
            x_values,
            y_values,
            color=COLOR_MAPPING[group],
            alpha=0.7,
            s=30,
            marker=MARKER_MAPPING[group],  # Use marker based on group
            linewidth=0,
        )

    # ------------ Wastewater plot styling---------- #

    ax_left.set_xscale("log")
    ax_left.set_xlabel("Relative abundance (wastewater)", fontsize=12)
    ax_left.set_yticks(panel_df["y"])
    ax_left.set_yticklabels(panel_df["species"], fontsize=10)
    ax_left.tick_params(axis="y", length=0)  # Remove y-axis tick marks
    ax_left.tick_params(axis="x", length=0)  # Remove x-axis tick marks
    ax_left.set_ylim(-0.5, max(panel_df["y"]) + 0.5)

    # Grid lines
    ax_left.grid(which="major", axis="x", linewidth=0.3, alpha=0.5)
    for y in panel_df["y"]:
        ax_left.axhline(
            y=y, color="lightgray", linestyle="-", linewidth=0.5, alpha=0.5, zorder=0
        )
    # X-limits
    max_ra = ww_data["relative_abundance"].max() * 2
    min_ra = 2e-10
    ax_left.set_xlim(min_ra, max_ra)

    # Major x-ticks
    major_ticks = [1e-9, 1e-8, 1e-7]
    zero_tick_location = ZERO_VALUE_REPLACEMENT
    all_ticks = [zero_tick_location] + major_ticks
    ax_left.set_xticks(all_ticks)
    labels = ["0"] + [f"$10^{{{int(np.log10(tick))}}}$" for tick in major_ticks]
    ax_left.set_xticklabels(labels)

    # Minor x-ticks
    minor_ticks = []
    for exp in range(-9, int(np.log10(max_ra)) + 1):
        base = 10**exp
        for mult in [2, 3, 4, 5, 6, 7, 8, 9]:
            tick_val = base * mult
            if min_ra <= tick_val <= max_ra:
                minor_ticks.append(tick_val)

    ax_left.set_xticks(minor_ticks, minor=True)

    # ---------- Spine styling ---------- #
    # Hide spines
    for spine in ["top", "right", "left", "bottom"]:
        ax_left.spines[spine].set_visible(False)

    # X-axis dashed between zero tick and 10^-9 to show discontinuity
    y_axis_bottom = ax_left.get_ylim()[0]
    min_nonzero_ra = ww_data[ww_data["relative_abundance"] > 0]["relative_abundance"].min()
    ax_left.plot(
        [min_ra, min_nonzero_ra],
        [y_axis_bottom, y_axis_bottom],
        color="black",
        linewidth=2,
        linestyle=(0, (1, 5)),
        solid_capstyle="butt",
    )

    # X-axis solid segment from smallest non-zero value to max_ra
    ax_left.plot(
        [min_nonzero_ra, max_ra],
        [y_axis_bottom, y_axis_bottom],
        color="black",
        linewidth=2,
        solid_capstyle="butt",
    )

    # -------------- Swab plot -------------- #
    bars = ax_right.barh(
        panel_df["y"],
        panel_df["n_positive_pools"],
        height=0.6,
        color=[COLOR_MAPPING[g] for g in panel_df["group"]],
    )

    # ----------- Swab plot styling --------- #
    ax_right.set_xlabel("Positive swab pools", fontsize=12)
    ax_right.set_xlim(0, max(panel_df["n_positive_pools"].max(), 1) * 1.1)
    ax_right.set_yticks([])  # Hide duplicate labels
    ax_right.tick_params(axis="x", length=0)  # Remove x-axis tick marks
    ax_right.set_ylim(ax_left.get_ylim())
    ax_right.grid(axis="x", linewidth=0.3, alpha=0.5, zorder=-5)

    # Set 0 values when no positive swab pools
    for bar in bars:
        width = bar.get_width()
        if width == 0:
            ax_right.text(
                width + 0.1,
                bar.get_y() + bar.get_height() / 2,
                f"{int(width)}",
                ha="left",
                va="center",
                fontsize=8,
            )

    for spine in ["top", "right"]:
        ax_right.spines[spine].set_visible(False)

    # -------------- Legend -------------- #
    present_groups = panel_df["group"].unique()
    ordered_groups = [group for group in GROUP_ORDER if group in present_groups]

    # Create legend entries in the specified order
    handles = []
    labels = []
    for group in ordered_groups:
        color = COLOR_MAPPING[group]
        marker = MARKER_MAPPING[group]
        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker=marker,
                color=color,
                markersize=8,
                linestyle="None",
                alpha=1.0,
            )
        )
        labels.append(group)

    # Position legend in the reserved bottom space
    fig.legend(
        handles,
        labels,
        loc="center",
        bbox_to_anchor=(0.4, 0.1),  # Position in the reserved bottom space
        ncol=len(ordered_groups) if len(ordered_groups) <= 4 else 3,
        frameon=False,
        fontsize=10,
    )

    # Save figure
    os.makedirs("figures", exist_ok=True)
    plt.savefig("figures/ra_and_presence.png", dpi=300)
    plt.savefig("figures/ra_and_presence.svg")


if __name__ == "__main__":
    main()
