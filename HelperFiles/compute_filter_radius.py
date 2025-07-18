#!/usr/bin/env python3
"""
Script to compute and visualize the distribution of edge lengths from a histogram CSV file,
and to suggest filter radii for topology optimization based on multiples of the mean edge length.

The script reads a CSV containing binned edge length data (start and end of bin, and count),
computes a weighted mean edge length, calculates suggested filter radii based on specified
multipliers, and plots the distribution with relevant annotations.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def compute_and_plot_filter_radius_from_csv(file_path, multiplier_values=[1.5, 2.0, 2.5, 3.0]):
    """
    Compute the mean edge length from a CSV containing histogram bins, and plot the histogram 
    with suggested filter radii marked.

    Parameters:
    -----------
    file_path : str
        Path to the CSV file containing the histogram of edge lengths.
        Expected format: index columns are bin start and end; the first column is the count.
    multiplier_values : list of float, optional
        Multipliers to apply to the mean edge length to suggest filter radii.
        Default is [1.5, 2.0, 2.5, 3.0].

    Returns:
    --------
    None
        Saves a PNG file of the plot and prints the computed values to the console.
    """

    # Load CSV and parse bin edges and counts
    df = pd.read_csv(file_path, header=0, index_col=[0, 1])

    bin_edges = []
    counts = []

    for idx, count in df.iloc[:, 0].items():
        try:
            left = float(idx[0])
            right = float(idx[1])
            bin_center = 0.5 * (left + right)
            bin_edges.append(bin_center)
            counts.append(count)
        except ValueError:
            continue  # Skip malformed rows

    bin_edges = np.array(bin_edges)
    counts = np.array(counts)

    # Compute weighted mean edge length
    mean_edge_length = np.average(bin_edges, weights=counts)

    # Compute filter radii as multiples of mean edge length
    filter_radii = {f"{m}×": m * mean_edge_length for m in multiplier_values}

    # Plotting
    plt.figure(figsize=(8, 5))
    plt.bar(bin_edges, counts, width=0.015, edgecolor='k', align='center', alpha=0.8)
    plt.axvline(x=mean_edge_length, color='red', linestyle='--', label=f'Mean = {mean_edge_length:.3f} m')

    for label, value in filter_radii.items():
        plt.axvline(x=value, linestyle=':', label=f'{label}Mean = {value:.3f} m')

    plt.xlabel("Edge Length (m)")
    plt.ylabel("Count")
    plt.title("Edge Length Distribution and Filter Radius Suggestions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("EdgeLengthDist.png", dpi=300)

    # Print results
    print(f"Mean edge length: {mean_edge_length:.5f} m")
    for label, value in filter_radii.items():
        print(f"Filter radius ({label}): {value:.5f} m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute filter radius from edge length histogram CSV")
    parser.add_argument("csv_path", type=str, help="Path to the edge length histogram CSV file")
    args = parser.parse_args()

    compute_and_plot_filter_radius_from_csv(args.csv_path)
