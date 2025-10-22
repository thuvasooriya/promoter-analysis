from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


class PromoterVisualizer:
    def __init__(self, results_dir="results"):
        self.results_dir = Path(results_dir)
        self.figures_dir = self.results_dir / "figures"
        self.figures_dir.mkdir(exist_ok=True)

    def load_data(self):
        self.ppm_df = pd.read_csv(
            self.results_dir / "task1_ppm" / "position_probability_matrix.csv",
            index_col=0,
        )
        self.results_df = pd.read_csv(
            self.results_dir / "task2_alignment" / "statistical_alignment_results.csv"
        )
        self.promoters_df = pd.read_csv(
            self.results_dir / "task1_ppm" / "manual_promoters.csv"
        )

    def plot_ppm_heatmap(self):
        fig, ax = plt.subplots(figsize=(8, 4))

        sns.heatmap(
            self.ppm_df.T,
            annot=True,
            fmt=".3f",
            cmap="YlOrRd",
            cbar_kws={"label": "Probability"},
            ax=ax,
        )

        ax.set_xlabel("Position in Promoter", fontsize=12)
        ax.set_ylabel("Nucleotide", fontsize=12)

        plt.tight_layout()
        plt.savefig(self.figures_dir / "ppm_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved: {self.figures_dir / 'ppm_heatmap.png'}")

    def plot_sequence_logo(self):
        fig, ax = plt.subplots(figsize=(10, 3))

        logo_df = self.ppm_df.copy()
        logo_df.index = range(len(logo_df))
        logo = logomaker.Logo(logo_df, ax=ax, color_scheme="classic")

        ax.set_xlabel("Position", fontsize=12)
        ax.set_ylabel("Probability", fontsize=12)
        # ax.set_title(
        #     "Sequence Logo - Consensus: TATAAT", fontsize=14, fontweight="bold"
        # )

        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "sequence_logo.png", dpi=300, bbox_inches="tight"
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'sequence_logo.png'}")

    def plot_position_probabilities(self):
        fig, ax = plt.subplots(figsize=(10, 4))

        positions = range(1, len(self.ppm_df) + 1)
        width = 0.8

        for base, color in [
            ("A", "#00CC00"),
            ("T", "#FF0000"),
            ("G", "#FFB300"),
            ("C", "#0000FF"),
        ]:
            heights = list(self.ppm_df[base].values)
            ax.bar(positions, heights, width, label=base, alpha=0.8, color=color)

        ax.set_xlabel("Position", fontsize=12)
        ax.set_ylabel("Probability", fontsize=12)
        ax.set_title(
            "Nucleotide Probability Distribution by Position\nConsensus: TATAAT",
            fontsize=14,
            fontweight="bold",
        )
        ax.set_xticks(positions)
        ax.set_xticklabels([f"Pos {i}" for i in positions])
        ax.legend()
        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "position_probabilities.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'position_probabilities.png'}")

    def plot_score_distribution(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        with_promoter = self.results_df[self.results_df["has_promoter"] == True][
            "best_score"
        ]
        without_promoter = self.results_df[self.results_df["has_promoter"] == False][
            "best_score"
        ]

        ax1.hist(
            with_promoter,
            bins=30,
            alpha=0.7,
            color="green",
            label="With Promoter",
            edgecolor="black",
        )
        ax1.hist(
            without_promoter,
            bins=30,
            alpha=0.7,
            color="red",
            label="Without Promoter",
            edgecolor="black",
        )
        ax1.axvline(
            with_promoter.min(),
            color="blue",
            linestyle="--",
            linewidth=2,
            label="Threshold",
        )
        ax1.set_xlabel("Alignment Score", fontsize=12)
        ax1.set_ylabel("Frequency", fontsize=12)
        ax1.set_title(
            "Score Distribution - Promoter vs Non-Promoter",
            fontsize=13,
            fontweight="bold",
        )
        ax1.legend()
        ax1.grid(axis="y", alpha=0.3)

        scores = self.results_df["best_score"]
        ax2.boxplot(
            [with_promoter, without_promoter],
            labels=["With Promoter", "Without Promoter"],
        )
        ax2.set_ylabel("Alignment Score", fontsize=12)
        ax2.set_title("Score Comparison - Box Plot", fontsize=13, fontweight="bold")
        ax2.grid(axis="y", alpha=0.3)

        plt.suptitle(
            f"Statistical Alignment Results (n={len(self.results_df)})",
            fontsize=14,
            fontweight="bold",
            y=1.02,
        )
        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "score_distributions.png", dpi=300, bbox_inches="tight"
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'score_distributions.png'}")

    def plot_detection_summary(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        detected = self.results_df["has_promoter"].sum()
        not_detected = len(self.results_df) - detected

        colors = ["#2ecc71", "#e74c3c"]
        ax1.pie(
            [detected, not_detected],
            labels=["Promoter Detected", "No Promoter"],
            autopct="%1.1f%%",
            startangle=90,
            colors=colors,
            explode=(0.05, 0),
        )
        ax1.set_title(
            "Promoter Detection Rate\n(1000 Test Sequences)",
            fontsize=13,
            fontweight="bold",
        )

        position_series = self.results_df[self.results_df["has_promoter"] == True][
            "best_position"
        ]
        position_counts = position_series.value_counts().sort_index()
        ax2.bar(
            list(position_counts.index),
            list(position_counts.values),
            color="#3498db",
            edgecolor="black",
        )
        ax2.set_xlabel("Position in Upstream Region", fontsize=12)
        ax2.set_ylabel("Count", fontsize=12)
        ax2.set_title(
            "Promoter Position Distribution\n(Where promoters were found)",
            fontsize=13,
            fontweight="bold",
        )
        ax2.grid(axis="y", alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "detection_summary.png", dpi=300, bbox_inches="tight"
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'detection_summary.png'}")

    def plot_training_data_analysis(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        w_contents = self.promoters_df["w_content"]
        w_mean = float(w_contents.mean())
        ax1.hist(w_contents, bins=15, color="#9b59b6", edgecolor="black", alpha=0.7)
        ax1.axvline(
            w_mean,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {w_mean:.2f}",
        )
        ax1.set_xlabel("W Content Score", fontsize=12)
        ax1.set_ylabel("Frequency", fontsize=12)
        ax1.set_title(
            "W Content Distribution in Training Set", fontsize=13, fontweight="bold"
        )
        ax1.legend()
        ax1.grid(axis="y", alpha=0.3)

        from collections import Counter

        sequences = list(self.promoters_df["promoter_sequence"])
        seq_counts = Counter(sequences)
        top_10 = dict(sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)[:10])

        ax2.barh(
            list(top_10.keys()),
            list(top_10.values()),
            color="#e67e22",
            edgecolor="black",
        )
        ax2.set_xlabel("Frequency", fontsize=12)
        ax2.set_ylabel("Promoter Sequence", fontsize=12)
        ax2.set_title(
            "Top 10 Most Common Promoter Sequences", fontsize=13, fontweight="bold"
        )
        ax2.grid(axis="x", alpha=0.3)

        # plt.suptitle(
        #     f"Training Data Analysis (n={len(self.promoters_df)})",
        #     fontsize=14,
        #     fontweight="bold",
        #     y=1.02,
        # )
        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "training_data_analysis.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'training_data_analysis.png'}")

    def plot_cross_validation_comparison(self):
        import json

        cv_file = (
            self.results_dir
            / "task3_cross_validation"
            / "cross_validation_results.json"
        )
        if not cv_file.exists():
            print(f"Cross-validation results not found: {cv_file}")
            return

        with open(cv_file) as f:
            cv_data = json.load(f)

        students = []
        detection_rates = []

        for student_id, results in cv_data.items():
            students.append(student_id)
            detection_rates.append(results["detection_rate"] * 100)

        task2_file = self.results_dir / "task2_alignment" / "task2_summary.json"
        if task2_file.exists():
            with open(task2_file) as f:
                task2_data = json.load(f)
            students.append("210657G\n(own)")
            detection_rates.append(task2_data["detection_rate"] * 100)

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = ["#3498db"] * (len(students) - 1) + ["#e74c3c"]
        bars = ax.bar(
            students, detection_rates, color=colors, edgecolor="black", alpha=0.8
        )

        mean_rate = sum(detection_rates) / len(detection_rates)
        ax.axhline(
            mean_rate,
            color="green",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {mean_rate:.2f}%",
        )

        ax.set_xlabel("Student Genome", fontsize=12, fontweight="bold")
        ax.set_ylabel("Detection Rate (%)", fontsize=12, fontweight="bold")
        # ax.set_title(
        #     "Cross-Validation: PPM Performance Across Genomes\nStudent 210657G's PPM",
        #     fontsize=14,
        #     fontweight="bold",
        # )
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        for i, (student, rate) in enumerate(zip(students, detection_rates)):
            ax.text(i, rate + 1, f"{rate:.1f}%", ha="center", va="bottom", fontsize=10)

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(
            self.figures_dir / "cross_validation_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        print(f"Saved: {self.figures_dir / 'cross_validation_comparison.png'}")

    def generate_all_plots(self):
        print("Generating visualizations...")
        self.load_data()
        self.plot_ppm_heatmap()
        self.plot_sequence_logo()
        self.plot_score_distribution()
        self.plot_detection_summary()
        self.plot_training_data_analysis()
        self.plot_cross_validation_comparison()
        print(f"\nAll visualizations saved to: {self.figures_dir}")


if __name__ == "__main__":
    visualizer = PromoterVisualizer()
    visualizer.generate_all_plots()
