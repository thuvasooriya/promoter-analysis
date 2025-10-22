"""
Run background distribution analysis
Compare training, test, and background distributions
Determine proper statistical threshold
"""

import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import sys
sys.path.append("src")

from background_distribution import BackgroundDistribution
from data_parser import GenomeParser
from ppm_builder import PPMBuilder

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("background_analysis.log"),
        logging.StreamHandler(),
    ],
)


def load_existing_results():
    """Load existing PPM and scores"""
    ppm_path = Path("results/task1_ppm/position_probability_matrix.csv")
    training_path = Path("results/task1_ppm/manual_promoters.csv")
    test_path = Path("results/task2_alignment/statistical_alignment_results.csv")

    if not all([ppm_path.exists(), training_path.exists(), test_path.exists()]):
        raise FileNotFoundError("Run tasks 1-2 first to generate PPM and scores")

    ppm_df = pd.read_csv(ppm_path, index_col=0)
    training_df = pd.read_csv(training_path)
    test_df = pd.read_csv(test_path)

    logging.info(f"Loaded PPM: {ppm_df.shape}")
    logging.info(f"Loaded training: {len(training_df)} sequences")
    logging.info(f"Loaded test: {len(test_df)} sequences")

    return ppm_df, training_df, test_df


def calculate_genome_gc_content(fasta_path: str) -> float:
    """Calculate genome GC content"""
    from Bio import SeqIO

    record = SeqIO.read(fasta_path, "fasta")
    seq = str(record.seq).upper()
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    logging.info(f"Genome GC content: {gc:.3f}")
    return gc


def run_background_analysis():
    """Run comprehensive background distribution analysis"""
    logging.info("=" * 60)
    logging.info("Background Distribution Analysis")
    logging.info("=" * 60)

    ppm_df, training_df, test_df = load_existing_results()

    data_dir = Path("data/210657G_GCA_900637025.1/GCA_900637025.1")
    fasta_path = str(data_dir / "GCA_900637025.1_46338_H01_genomic.fna")
    gff_path = str(data_dir / "genomic.gff")

    genome_gc = calculate_genome_gc_content(fasta_path)

    bg_analyzer = BackgroundDistribution(ppm_df)

    # Method 1: Random sequences with genome GC content
    logging.info("\n" + "=" * 60)
    logging.info("Method 1: Random Sequences (Genome GC Content)")
    logging.info("=" * 60)
    bg_random_genome = bg_analyzer.analyze_background(
        method="random", n_sequences=10000, gc_content=genome_gc
    )

    # Method 2: Random sequences with uniform distribution
    logging.info("\n" + "=" * 60)
    logging.info("Method 2: Random Sequences (Uniform GC=0.5)")
    logging.info("=" * 60)
    bg_random_uniform = bg_analyzer.analyze_background(
        method="random", n_sequences=10000, gc_content=0.5
    )

    # Method 3: Intergenic sequences
    logging.info("\n" + "=" * 60)
    logging.info("Method 3: Intergenic Sequences (Real Genome)")
    logging.info("=" * 60)
    bg_intergenic = bg_analyzer.analyze_background(
        method="intergenic",
        n_sequences=10000,
        fasta_path=fasta_path,
        gff_path=gff_path,
    )

    # Re-score training promoters using PPM (CSV contains W-content, not PPM scores)
    from scoring import SequenceScorer
    
    scorer = SequenceScorer(ppm_df)
    training_scores = [
        scorer.score_sequence(seq, normalize=True)
        for seq in training_df["promoter_sequence"]
    ]
    test_scores = test_df["best_score"].tolist()
    
    logging.info(f"Re-scored {len(training_scores)} training promoters with PPM")

    # Compare distributions
    logging.info("\n" + "=" * 60)
    logging.info("Distribution Comparison")
    logging.info("=" * 60)

    comparison = bg_analyzer.compare_distributions(
        training_scores, test_scores, bg_intergenic["scores"]
    )

    # Calculate recommended thresholds
    logging.info("\n" + "=" * 60)
    logging.info("Recommended Thresholds")
    logging.info("=" * 60)

    thresholds = {
        "training_mean_minus_2std": np.mean(training_scores)
        - 2 * np.std(training_scores),
        "background_mean_plus_2std": bg_intergenic["threshold_2std"],
        "background_mean_plus_3std": bg_intergenic["threshold_3std"],
        "background_95pct": bg_intergenic["threshold_95pct"],
        "background_99pct": bg_intergenic["threshold_99pct"],
        "static_current": -10.0,
    }

    logging.info("Threshold Comparison:")
    for name, value in thresholds.items():
        count = sum(1 for s in test_scores if s > value)
        rate = count / len(test_scores) * 100
        logging.info(f"  {name:30s}: {value:7.3f}  →  {count:4d}/1000 ({rate:5.1f}%)")

    # Save results
    results_dir = Path("results/background_analysis")
    results_dir.mkdir(exist_ok=True, parents=True)

    summary = {
        "genome_gc": genome_gc,
        "training_stats": comparison["training"],
        "test_stats": comparison["test"],
        "background_methods": {
            "random_genome_gc": {
                "gc_content": bg_random_genome["gc_content"],
                "mean": bg_random_genome["mean"],
                "std": bg_random_genome["std"],
                "threshold_2std": bg_random_genome["threshold_2std"],
                "threshold_3std": bg_random_genome["threshold_3std"],
            },
            "random_uniform": {
                "gc_content": bg_random_uniform["gc_content"],
                "mean": bg_random_uniform["mean"],
                "std": bg_random_uniform["std"],
                "threshold_2std": bg_random_uniform["threshold_2std"],
                "threshold_3std": bg_random_uniform["threshold_3std"],
            },
            "intergenic": {
                "gc_content": bg_intergenic["gc_content"],
                "mean": bg_intergenic["mean"],
                "std": bg_intergenic["std"],
                "threshold_2std": bg_intergenic["threshold_2std"],
                "threshold_3std": bg_intergenic["threshold_3std"],
            },
        },
        "thresholds": thresholds,
        "test_detection_rates": {
            name: sum(1 for s in test_scores if s > value) / len(test_scores)
            for name, value in thresholds.items()
        },
    }

    with open(results_dir / "background_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Create visualizations
    create_visualizations(
        training_scores,
        test_scores,
        bg_intergenic["scores"],
        bg_random_genome["scores"],
        thresholds,
        results_dir,
    )

    logging.info(f"\nResults saved to {results_dir}")

    return summary


def create_visualizations(
    training_scores, test_scores, bg_intergenic, bg_random, thresholds, output_dir
):
    """Create comparison visualizations"""
    logging.info("\nGenerating visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: All distributions
    ax = axes[0, 0]
    ax.hist(
        training_scores,
        bins=30,
        alpha=0.7,
        label=f"Training (n={len(training_scores)})",
        color="green",
    )
    ax.hist(
        test_scores, bins=50, alpha=0.5, label=f"Test (n={len(test_scores)})", color="blue"
    )
    ax.hist(
        bg_intergenic,
        bins=50,
        alpha=0.5,
        label=f"Background Intergenic (n={len(bg_intergenic)})",
        color="red",
    )
    ax.axvline(thresholds["static_current"], color="black", linestyle="--", label="Static (-10.0)")
    ax.axvline(
        thresholds["background_mean_plus_2std"],
        color="orange",
        linestyle="--",
        label="BG μ+2σ",
    )
    ax.set_xlabel("Score")
    ax.set_ylabel("Frequency")
    ax.set_title("Distribution Comparison: Training vs Test vs Background")
    ax.legend()
    ax.grid(alpha=0.3)

    # Plot 2: Test vs Background only (zoomed)
    ax = axes[0, 1]
    ax.hist(
        test_scores, bins=50, alpha=0.7, label=f"Test (n={len(test_scores)})", color="blue"
    )
    ax.hist(
        bg_intergenic,
        bins=50,
        alpha=0.5,
        label=f"Background (n={len(bg_intergenic)})",
        color="red",
    )
    ax.axvline(thresholds["static_current"], color="black", linestyle="--", label="Static (-10.0)")
    ax.axvline(
        thresholds["background_mean_plus_2std"],
        color="orange",
        linestyle="--",
        label="BG μ+2σ",
    )
    ax.set_xlabel("Score")
    ax.set_ylabel("Frequency")
    ax.set_title("Test vs Background (Overlapping Ranges)")
    ax.legend()
    ax.grid(alpha=0.3)

    # Plot 3: Threshold comparison
    ax = axes[1, 0]
    threshold_names = list(thresholds.keys())
    threshold_values = list(thresholds.values())
    detection_rates = [
        sum(1 for s in test_scores if s > t) / len(test_scores) * 100
        for t in threshold_values
    ]

    colors = ["gray" if "static" in name else "blue" for name in threshold_names]
    bars = ax.barh(threshold_names, detection_rates, color=colors, alpha=0.7)
    ax.axvline(33.6, color="red", linestyle="--", label="Current (33.6%)")
    ax.set_xlabel("Detection Rate (%)")
    ax.set_title("Detection Rates by Threshold Method")
    ax.legend()
    ax.grid(alpha=0.3, axis="x")

    # Plot 4: Box plots
    ax = axes[1, 1]
    data_to_plot = [training_scores, test_scores, bg_intergenic, bg_random]
    labels = ["Training", "Test", "BG Intergenic", "BG Random"]
    bp = ax.boxplot(
        data_to_plot, labels=labels, patch_artist=True, showmeans=True, widths=0.6
    )

    colors = ["green", "blue", "red", "orange"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.axhline(thresholds["static_current"], color="black", linestyle="--", alpha=0.5)
    ax.set_ylabel("Score")
    ax.set_title("Score Distribution Comparison (Box Plots)")
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_dir / "background_analysis.png", dpi=300, bbox_inches="tight")
    logging.info(f"Saved visualization: {output_dir / 'background_analysis.png'}")

    # Additional plot: ROC-style curve
    fig, ax = plt.subplots(figsize=(10, 8))

    threshold_range = np.linspace(-30, -5, 100)
    detection_rates_curve = [
        sum(1 for s in test_scores if s > t) / len(test_scores) * 100
        for t in threshold_range
    ]

    ax.plot(threshold_range, detection_rates_curve, linewidth=2, color="blue")

    for name, value in thresholds.items():
        rate = sum(1 for s in test_scores if s > value) / len(test_scores) * 100
        if "static" in name:
            ax.plot(value, rate, "ro", markersize=10, label=f"{name} ({rate:.1f}%)")
        elif "background" in name:
            ax.plot(value, rate, "go", markersize=8, alpha=0.7, label=f"{name} ({rate:.1f}%)")

    ax.set_xlabel("Threshold", fontsize=12)
    ax.set_ylabel("Detection Rate (%)", fontsize=12)
    ax.set_title("Detection Rate vs Threshold", fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "threshold_curve.png", dpi=300, bbox_inches="tight")
    logging.info(f"Saved visualization: {output_dir / 'threshold_curve.png'}")

    plt.close("all")


if __name__ == "__main__":
    results = run_background_analysis()
    logging.info("\n✓ Background analysis complete")
