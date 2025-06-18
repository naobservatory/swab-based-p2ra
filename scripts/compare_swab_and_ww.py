#!/usr/bin/env python3

import csv
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from datetime import datetime
import matplotlib.pyplot as plt

GENOME_LENGTH = 7000
DOUBLING_TIME = 4
SHEDDING_DURATION = 7
POOL_SIZE = 500
HV_COVERAGE_SWAB = 2000
HV_COVERAGE_WW = 300
SEQ_DEPTH_SWAB = 243294
SEQ_DEPTH_WW = 3000000000

results_swab = []
results_ww = []

posterior_df = pd.read_csv("tables/posteriors.tsv", sep="\t")
target_viruses = ["Rhinovirus A", "Rhinovirus B", "Rhinovirus C"]


fig, axs = plt.subplots(len(target_viruses), 1, figsize=(10, 10))

for plot_index, virus in enumerate(target_viruses):
    rhino_df = posterior_df[posterior_df["species"] == virus]
    mu_swabs = rhino_df["mu_swab"].values
    mu_wws = rhino_df["mu_ww"].values

    for i in range(5000):
        incidence = 1e-6  # 1 in 1 million
        cumulative_incidence = 0
        incidences = []

        swab_detected = False
        ww_detected = False

        for sampling_day in range(60):  # We sample every day, for 60 days
            if swab_detected and ww_detected:
                break

            incidence = incidence * 2 ** (1 / DOUBLING_TIME)
            cumulative_incidence += incidence
            if incidence > 1:
                break

            incidences.append(incidence)
            prevalence = sum(
                incidences[-SHEDDING_DURATION:]
            )  # Given shedding duration, each sample can contain samples from people infected in the last 7 days

            if not swab_detected:
                n_positive = np.random.binomial(POOL_SIZE, prevalence)
                ra_swabs = []
                if n_positive == 0:
                    ra_swabs = []
                if n_positive == 1:
                    ra_swab = np.random.choice(mu_swabs)
                    ra_swabs = [ra_swab]
                else:
                    for _ in range(n_positive):
                        ra_swab = np.random.choice(mu_swabs)
                        ra_swabs.append(ra_swab)

                summed_ra_swabs = sum(ra_swabs)

                sample_ra = summed_ra_swabs / POOL_SIZE

                pos_reads = round(sample_ra * SEQ_DEPTH_SWAB)

                for read in range(pos_reads):
                    detection_likelihood = HV_COVERAGE_SWAB / GENOME_LENGTH
                    detected = np.random.random() < detection_likelihood
                    if detected:
                        results_swab.append(cumulative_incidence)
                        swab_detected = True
                        break

            if not ww_detected:
                ra_ww = np.random.choice(mu_wws)
                ra_ww = ra_ww / 100  # Converting to RA(1%)
                ra_ww = ra_ww * prevalence
                pos_reads = round(ra_ww * SEQ_DEPTH_WW)
                detection_likelihood = HV_COVERAGE_WW / GENOME_LENGTH
                for read in range(pos_reads):
                    detected = np.random.random() < detection_likelihood
                    if detected:
                        results_ww.append(cumulative_incidence)
                        ww_detected = True
                        break

    results_sorted_swab = sorted(results_swab)
    x_percentages_swab = [
        i / len(results_sorted_swab) * 100 for i in range(len(results_sorted_swab))
    ]

    results_sorted_ww = sorted(results_ww)
    x_percentages_ww = [
        i / len(results_sorted_ww) * 100 for i in range(len(results_sorted_ww))
    ]

    axs[plot_index].plot(x_percentages_swab, results_sorted_swab)
    axs[plot_index].plot(x_percentages_ww, results_sorted_ww)
    axs[plot_index].set_ylabel("Cumulative Incidence (%)")
    axs[plot_index].yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, p: f"{x*100:.1f}%")
    )

    axs[plot_index].set_ylim(0, 0.1)

    if plot_index == len(target_viruses) - 1:
        axs[plot_index].set_xlabel("Percentile")
    axs[plot_index].set_title(virus)
    axs[plot_index].legend(["Swab", "WW"])
    axs[plot_index].grid(True, alpha=0.3)
plt.show()
