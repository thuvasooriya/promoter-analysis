"""
Fixed statistical alignment implementation
"""

import logging
from typing import Dict, List

import numpy as np
import pandas as pd


class StatisticalAligner:
    def __init__(self, ppm_df: pd.DataFrame):
        """Initialize with Position Probability Matrix"""
        self.ppm_df = ppm_df
        self.ppm_length = len(ppm_df)
        self.consensus_score = self._calculate_consensus_score()
        self.threshold = -10.0
        logging.info(
            f"Initialized StatisticalAligner with PPM of length {self.ppm_length}"
        )

    def _calculate_consensus_score(self) -> float:
        """Calculate consensus sequence benchmark score"""
        consensus_probs = []
        for _, row in self.ppm_df.iterrows():
            max_prob = row.max()
            consensus_probs.append(np.log(max_prob))
        return sum(consensus_probs)

    def score_sequence(self, sequence: str) -> float:
        """Score a sequence using PPM"""
        if len(sequence) != self.ppm_length:
            raise ValueError(f"Sequence length must be {self.ppm_length}")

        log_score = 0.0
        for pos, base in enumerate(sequence):
            if base in ["A", "C", "G", "T"]:
                prob = self.ppm_df.iloc[pos][base]
                if prob > 0:
                    log_score += np.log(prob)
                else:
                    log_score += np.log(1e-10)  # Small probability for zero values

        return log_score

    def sliding_window_analysis(self, sequence: str) -> List[Dict]:
        """Perform sliding window analysis on longer sequence"""
        results = []

        if len(sequence) < self.ppm_length:
            return results

        for i in range(len(sequence) - self.ppm_length + 1):
            subseq = sequence[i : i + self.ppm_length]
            score = self.score_sequence(subseq)
            normalized_score = score - self.consensus_score

            results.append({"position": i, "sequence": subseq, "score": normalized_score})

        return results

    def set_threshold(self, training_scores: List[float], method="mean_minus_2std"):
        """Set detection threshold based on training data"""
        if method == "mean_minus_2std":
            mean_score = np.mean(training_scores)
            std_score = np.std(training_scores)
            self.threshold = mean_score - 2 * std_score
        elif method == "percentile_5":
            self.threshold = np.percentile(training_scores, 5)
        else:
            self.threshold = np.mean(training_scores) - np.std(training_scores)

        logging.info(f"Set detection threshold to {self.threshold:.3f}")

    def analyze_upstream_regions(self, upstream_regions: List[Dict]) -> List[Dict]:
        """Analyze upstream regions for promoters"""
        results = []

        logging.info(f"Analyzing {len(upstream_regions)} upstream regions")

        for region in upstream_regions:
            sequence = region["upstream_sequence"]

            if len(sequence) < self.ppm_length:
                continue

            if len(sequence) == self.ppm_length:
                score = self.score_sequence(sequence)
                normalized_score = score - self.consensus_score

                results.append(
                    {
                        "gene_id": region["gene_id"],
                        "upstream_sequence": sequence,
                        "best_score": normalized_score,
                        "best_position": 0,
                        "best_subsequence": sequence,
                        "has_promoter": normalized_score > self.threshold,
                    }
                )
            else:
                window_results = self.sliding_window_analysis(sequence)

                if window_results:
                    best_result = max(window_results, key=lambda x: x["score"])

                    results.append(
                        {
                            "gene_id": region["gene_id"],
                            "upstream_sequence": sequence,
                            "best_score": best_result["score"],
                            "best_position": best_result["position"],
                            "best_subsequence": best_result["sequence"],
                            "has_promoter": best_result["score"] > self.threshold,
                        }
                    )

        promoter_count = sum(1 for r in results if r["has_promoter"])
        logging.info(
            f"Found {promoter_count}/{len(results)} sequences with predicted promoters"
        )

        return results
