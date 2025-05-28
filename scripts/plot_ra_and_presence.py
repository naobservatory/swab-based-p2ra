#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import textwrap
from metadata_utils import first_level_mapping, second_level_mapping

# Define color palettes for different virus groups
color_mapping = {
    "Coronaviruses (seasonal)": "#05a4a5",
    "Coronaviruses (SARS-CoV-2)": "#445681",
    "Rhinoviruses": "#ba5c97",
    "Mononegavirales": "#8CCEA4",
    "Influenza": "#E08F60",
    "Other": "#9D9D9D"
}

markers = ['s', 'o', '^']  # square, triangle, circle

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
    grid = df[["date", "location", "all_reads"]].drop_duplicates().merge(
        df[["species", "group"]].drop_duplicates(), how="cross"
    )

    # Merge with actual data, filling zeros where no detection occurred
    df_complete = (
        grid.merge(df, how="left")
        .fillna({"dedup_hv": 0, "non_dedup_hv": 0})
    )

    # Calculate relative abundance using dedup_hv / all_reads
    df_complete["relative_abundance"] = df_complete["dedup_hv"] / df_complete["all_reads"]

    return df_complete

# Load swab pathogen presence data
def load_swab_presence_data():
    # Load swab data
    swabs = pd.read_csv("tables/swabs-ra-summary.tsv", sep="\t")
    # Drop assignment column and collapse by summing dedup_hv and all_reads
    swabs_df = swabs.drop(columns="assignment").groupby(['date', 'location', 'species', 'group']).agg({
        'dedup_hv': 'sum',
        'non_dedup_hv': 'sum',
        'all_reads': 'first'
    }).reset_index()

    # Create a presence indicator (1 if virus is present)
    swabs_df["present"] = (swabs_df["dedup_hv"] > 0).astype(int)

    return swabs_df



def handle_zero_values(values):
    values = values.to_numpy(float)
    mask = values == 0
    values[mask] = ZERO_VALUE_REPLACEMENT
    return values, mask

def apply_jitter_to_zeros(x, y, mask, x_scale=0.1, y_scale=0.2):
    if mask.any():
        factors = np.random.uniform(1 - x_scale, 1 + x_scale, mask.sum())
        x[mask] *= factors
        y[mask] += np.random.uniform(-y_scale, y_scale, mask.sum())
    return x, y


# ---------- Load & prepare data ---------- #
def main():
    ww = load_ww_abundance_data()
    sw = load_swab_presence_data()

    # Get all unique species from wastewater data
    all_species = ww[['species', 'group']].drop_duplicates()

    # Count positive swab pools
    swab_counts = (
        sw[sw['present'] == 1]
        .groupby(['species', 'group'])
        .size()
        .reset_index(name='n_positive_pools')
    )

    # Calculate median relative abundance per species
    median_ra = (
        ww.groupby('species')['relative_abundance']
        .median()
        .rename('median_ra')
        .reset_index()
    )

    # Calculate group average relative abundance
    group_avg_ra = (
        ww.groupby('group')['relative_abundance']
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
    print(panel_df)

    # Sort by group average RA (ascending = lowest abundance groups at top)
    # and within each group sort by species median RA (descending)
    panel_df = panel_df.sort_values(['group_avg_ra', 'median_ra'],
                                ascending=[False, False]
                                ).reset_index(drop=True)


    # Create a y-position lookup with explicit positions (ensures alignment)
    num_species = len(panel_df)
    spacing = 1.0  # Consistent spacing
    panel_df['y'] = np.arange(num_species) * spacing
    y_lookup = panel_df.set_index('species')['y'].to_dict()


    fig = plt.figure(figsize=(14, 7.2), dpi=450)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 0.15], width_ratios=[3, 1], wspace=0.05, hspace=0)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])

    # Define marker cycling for different groups
    group_to_marker = {}

    # Assign markers to groups
    for i, group in enumerate(panel_df['group'].unique()):
        group_to_marker[group] = markers[i % len(markers)]

    # Wastewater plot
    for _, row in panel_df.iterrows():
        sp = row['species']
        grp = row['group']
        y = row['y']
        sp_data = ww[ww['species'] == sp]
        # Handle zero values and apply jitter
        x_values, zero_mask = handle_zero_values(sp_data['relative_abundance'])
        y_values = np.full(len(sp_data), y)
        x_jittered, y_jittered = apply_jitter_to_zeros(x_values, y_values, zero_mask)

        ax_left.scatter(
            x_jittered,
            y_jittered,
            color=color_mapping.get(grp, 'gray'),
            alpha=0.7,  # No transparency
            s=30,
            marker=group_to_marker[grp],  # Use marker based on group
            linewidth=0
        )

    ax_left.set_xscale('log')
    ax_left.set_xlabel('Relative abundance (wastewater)', fontsize=12)
    # ax_left.set_title('Wastewater Abundance', fontsize=14)
    ax_left.set_yticks(panel_df['y'])
    ax_left.set_yticklabels(panel_df['species'], fontsize=10)
    ax_left.tick_params(axis='y', length=0)  # Remove y-axis tick marks
    ax_left.tick_params(axis='x', length=0)  # Remove x-axis tick marks
    ax_left.grid(False)  # Remove default grid
    # Set y-limits with explicit padding to match right axis
    ax_left.set_ylim(-0.5, max(panel_df['y']) + 0.5)

    # Only show grid lines at order of magnitude intervals
    ax_left.grid(which='major', axis='x', linewidth=0.3, alpha=0.5)
    ax_left.grid(which='minor', axis='x', linestyle='', alpha=0)  # Remove minor grid lines

    # Set x-limits
    max_ra = ww['relative_abundance'].max() * 2
    min_ra = 2e-10
    ax_left.set_xlim(min_ra, max_ra)

    # Set up custom ticks and labels
    from matplotlib.ticker import LogLocator, LogFormatter

    # Create custom tick locations
    major_ticks = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0]
    zero_tick_location = np.sqrt(1e-9 * 1e-10)  # Geometric mean ≈ 3.16e-10

    # Filter major ticks to only include those within our range
    major_ticks = [tick for tick in major_ticks if tick >= min_ra and tick <= max_ra]
    all_ticks = [zero_tick_location] + major_ticks

    # Set custom ticks
    ax_left.set_xticks(all_ticks)

    # Create custom labels
    labels = ['0'] + [f'$10^{{{int(np.log10(tick))}}}$' for tick in major_ticks]
    ax_left.set_xticklabels(labels)

    # Add minor ticks from 1e-9 to max, but not below 1e-9
    # Create manual minor tick locations
    minor_ticks = []
    for exp in range(-9, int(np.log10(max_ra)) + 1):  # From 1e-9 to max_ra
        base = 10**exp
        for mult in [2, 3, 4, 5, 6, 7, 8, 9]:
            tick_val = base * mult
            if min_ra <= tick_val <= max_ra and tick_val >= 1e-9:
                minor_ticks.append(tick_val)

    ax_left.set_xticks(minor_ticks, minor=True)
    ax_left.tick_params(which='minor', bottom=True, length=3)

    # Make the x-axis dashed between zero tick and 10^-9 to show discontinuity
    zero_tick_location = np.sqrt(1e-9 * 1e-10)  # ≈ 3.16e-10
    y_axis_bottom = ax_left.get_ylim()[0]

    # Hide the default bottom spine
    ax_left.spines['bottom'].set_visible(False)

    # Find the smallest non-zero relative abundance value
    min_nonzero_ra = ww[ww['relative_abundance'] > 0]['relative_abundance'].min()

    # Dashed segment from zero_tick_location to smallest non-zero value
    ax_left.plot([min_ra, min_nonzero_ra], [y_axis_bottom, y_axis_bottom],
                color='black', linewidth=2, linestyle=(0, (1, 5)), solid_capstyle='butt')

    # Solid segment from smallest non-zero value to max_ra
    ax_left.plot([min_nonzero_ra, max_ra], [y_axis_bottom, y_axis_bottom],
                color='black', linewidth=2, solid_capstyle='butt')

    # Add horizontal grid lines aligned with each virus
    for y in panel_df['y']:
        ax_left.axhline(y=y, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.5, zorder=0)

    # ---- Right: swab positive-count bars ---- #
    bars = ax_right.barh(
        panel_df['y'],
        panel_df['n_positive_pools'],
        height=0.6,
        color=[color_mapping.get(g, 'gray') for g in panel_df['group']]
    )

    ax_right.set_xlabel('Positive swab pools', fontsize=12)
    # ax_right.set_title('Swab Positivity', fontsize=14)
    ax_right.set_xlim(0, max(panel_df['n_positive_pools'].max(), 1) * 1.1)
    ax_right.set_yticks([])  # hide duplicate labels
    ax_right.tick_params(axis='x', length=0)  # Remove x-axis tick marks
    # Set y-limits exactly the same as left axis for perfect alignment
    ax_right.set_ylim(ax_left.get_ylim())
    ax_right.grid(axis='x', linewidth=0.3, alpha=0.5,zorder=-5)

    # Add value labels to bars for easier reading
    for bar in bars:
        width = bar.get_width()
        if width == 0:
            ax_right.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                        f'{int(width)}', ha='left', va='center', fontsize=8)

    # ---- Add legend for virus groups ---- #
    # Define the desired order explicitly based on what we can see in the data
    desired_order = [
        "Rhinoviruses",
        "Influenza",
        "Mononegavirales",


        "Coronaviruses (SARS-CoV-2)",
        "Coronaviruses (seasonal)",
    ]

    # Filter to include only groups that are actually in our data
    present_groups = panel_df['group'].unique()
    ordered_groups = [group for group in desired_order if group in present_groups]

    # Create legend entries in the specified order
    handles = []
    labels = []
    for group in ordered_groups:
        color = color_mapping.get(group, 'gray')
        marker = group_to_marker.get(group, 'o')
        handles.append(plt.Line2D([0], [0], marker=marker, color=color, markersize=8,
                            linestyle='None', alpha=1.0))
        labels.append(group)

    # ---- Final touches ---- #
    for spine in ['top', 'right', 'left']:  # Also remove the y-axis spine on the left plot (bottom already handled above)
        ax_left.spines[spine].set_visible(False)

    for spine in ['top', 'right']:
        ax_right.spines[spine].set_visible(False)

    # Position legend in the reserved bottom space
    fig.legend(handles, labels,
            loc='center',
            bbox_to_anchor=(0.4, 0.1),  # Position in the reserved bottom space
            ncol=len(ordered_groups) if len(ordered_groups) <= 4 else 3,
            frameon=False,
            fontsize=10)

    # Save figure
    plt.savefig('figures/wastewater_vs_swab.png')
    print("Saved to figures/wastewater_vs_swab.png")

if __name__ == "__main__":
    main()