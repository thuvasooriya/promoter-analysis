"""
Background distribution analysis for promoter detection threshold
Generates random genomic sequences to model non-promoter baseline
"""

import logging
import random
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

from scoring import SequenceScorer


class BackgroundDistribution:
    def __init__(self, ppm_df: pd.DataFrame):
        """Initialize with PPM for scoring"""
        self.ppm_df = ppm_df
        self.ppm_length = len(ppm_df)
        self.bases = ["A", "C", "G", "T"]
        self.scorer = SequenceScorer(ppm_df)
        self.consensus_score = self.scorer.consensus_score
        logging.info(f"Initialized background distribution analyzer")

    def score_sequence(self, sequence: str) -> float:
        """Score a sequence using PPM (normalized)"""
        return self.scorer.score_sequence(sequence, normalize=True)

    def generate_random_sequences(
        self, n_sequences: int, gc_content: float = 0.5
    ) -> List[str]:
        """Generate random sequences with specified GC content"""
        at_prob = (1 - gc_content) / 2
        gc_prob = gc_content / 2

        base_probs = {
            "A": at_prob,
            "T": at_prob,
            "C": gc_prob,
            "G": gc_prob,
        }

        sequences = []
        for _ in range(n_sequences):
            seq = "".join(
                random.choices(
                    self.bases, weights=[base_probs[b] for b in self.bases], k=self.ppm_length
                )
            )
            sequences.append(seq)

        logging.info(
            f"Generated {n_sequences} random {self.ppm_length}bp sequences (GC={gc_content:.2f})"
        )
        return sequences

    def extract_intergenic_sequences(
        self, fasta_path: str, gff_path: str, n_sequences: int = 10000
    ) -> List[str]:
        """Extract random sequences from intergenic regions"""
        # Parse genome
        genome_record = SeqIO.read(fasta_path, "fasta")
        genome_seq = str(genome_record.seq)
        genome_length = len(genome_seq)

        # Parse GFF to identify genic regions
        genic_regions = []
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) >= 9 and fields[2] == "CDS":
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    genic_regions.append((start, end))

        genic_regions.sort()
        logging.info(f"Found {len(genic_regions)} genic regions")

        # Extract random intergenic sequences
        sequences = []
        attempts = 0
        max_attempts = n_sequences * 10

        while len(sequences) < n_sequences and attempts < max_attempts:
            pos = random.randint(0, genome_length - self.ppm_length)
            seq = genome_seq[pos : pos + self.ppm_length].upper()

            if "N" in seq:
                attempts += 1
                continue

            # Check if overlaps with genic region (with buffer)
            buffer = 50
            in_genic = False
            for start, end in genic_regions:
                if pos + self.ppm_length + buffer >= start and pos - buffer <= end:
                    in_genic = True
                    break

            if not in_genic:
                sequences.append(seq)

            attempts += 1

        logging.info(
            f"Extracted {len(sequences)} intergenic sequences ({attempts} attempts)"
        )
        return sequences

    def analyze_background(
        self,
        method: str = "random",
        n_sequences: int = 10000,
        gc_content: float = 0.5,
        fasta_path: str = "",
        gff_path: str = "",
    ) -> Dict:
        """Analyze background distribution and calculate threshold"""
        if method == "random":
            sequences = self.generate_random_sequences(n_sequences, gc_content)
            actual_gc = gc_content
        elif method == "intergenic":
            if not fasta_path or not gff_path:
                raise ValueError("fasta_path and gff_path required for intergenic method")
            sequences = self.extract_intergenic_sequences(
                fasta_path, gff_path, n_sequences
            )
            actual_gc = self._calculate_gc(sequences)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Score all sequences
        scores = [self.score_sequence(seq) for seq in sequences]

        # Calculate statistics
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        median_score = np.median(scores)

        # Calculate thresholds at different confidence levels
        threshold_1std = mean_score + std_score
        threshold_2std = mean_score + 2 * std_score
        threshold_3std = mean_score + 3 * std_score
        threshold_95pct = np.percentile(scores, 95)
        threshold_99pct = np.percentile(scores, 99)

        results = {
            "method": method,
            "n_sequences": len(sequences),
            "gc_content": actual_gc,
            "mean": mean_score,
            "std": std_score,
            "median": median_score,
            "min": np.min(scores),
            "max": np.max(scores),
            "threshold_1std": threshold_1std,
            "threshold_2std": threshold_2std,
            "threshold_3std": threshold_3std,
            "threshold_95pct": threshold_95pct,
            "threshold_99pct": threshold_99pct,
            "scores": scores,
        }

        logging.info(f"\nBackground Distribution Analysis ({method}):")
        logging.info(f"  Mean: {mean_score:.3f}")
        logging.info(f"  Std:  {std_score:.3f}")
        logging.info(f"  Range: [{np.min(scores):.3f}, {np.max(scores):.3f}]")
        logging.info(f"  Threshold (μ + 1σ): {threshold_1std:.3f}")
        logging.info(f"  Threshold (μ + 2σ): {threshold_2std:.3f}")
        logging.info(f"  Threshold (μ + 3σ): {threshold_3std:.3f}")
        logging.info(f"  Threshold (95th %ile): {threshold_95pct:.3f}")
        logging.info(f"  Threshold (99th %ile): {threshold_99pct:.3f}")

        return results

    def _calculate_gc(self, sequences: List[str]) -> float:
        """Calculate GC content of sequences"""
        total_bases = 0
        gc_bases = 0
        for seq in sequences:
            total_bases += len(seq)
            gc_bases += seq.count("G") + seq.count("C")
        return gc_bases / total_bases if total_bases > 0 else 0.0

    def compare_distributions(
        self, training_scores: List[float], test_scores: List[float], background_scores: List[float]
    ) -> Dict:
        """Compare training, test, and background distributions"""
        comparison = {
            "training": {
                "mean": np.mean(training_scores),
                "std": np.std(training_scores),
                "median": np.median(training_scores),
                "min": np.min(training_scores),
                "max": np.max(training_scores),
                "n": len(training_scores),
            },
            "test": {
                "mean": np.mean(test_scores),
                "std": np.std(test_scores),
                "median": np.median(test_scores),
                "min": np.min(test_scores),
                "max": np.max(test_scores),
                "n": len(test_scores),
            },
            "background": {
                "mean": np.mean(background_scores),
                "std": np.std(background_scores),
                "median": np.median(background_scores),
                "min": np.min(background_scores),
                "max": np.max(background_scores),
                "n": len(background_scores),
            },
        }

        logging.info("\nDistribution Comparison:")
        for dist_name, stats in comparison.items():
            logging.info(
                f"  {dist_name.capitalize()}: μ={stats['mean']:.2f}, σ={stats['std']:.2f}, "
                f"range=[{stats['min']:.2f}, {stats['max']:.2f}], n={stats['n']}"
            )

        return comparison
