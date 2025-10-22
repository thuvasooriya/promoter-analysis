"""
Fixed statistical alignment implementation
"""

import logging
from typing import Dict, List

import numpy as np
import pandas as pd

from scoring import SequenceScorer


class StatisticalAligner:
    def __init__(self, ppm_df: pd.DataFrame, threshold_method: str = "empirical"):
        """Initialize with Position Probability Matrix"""
        self.ppm_df = ppm_df
        self.ppm_length = len(ppm_df)
        self.scorer = SequenceScorer(ppm_df)
        self.consensus_score = self.scorer.consensus_score
        self.threshold_method = threshold_method
        
        if threshold_method == "consensus":
            self.threshold = 0.0
        else:
            self.threshold = -10.0
        
        logging.info(
            f"Initialized StatisticalAligner with PPM of length {self.ppm_length}"
        )
        logging.info(f"Threshold method: {threshold_method}, threshold: {self.threshold}")

    def score_sequence(self, sequence: str) -> float:
        """Score a sequence using PPM (normalized)"""
        return self.scorer.score_sequence(sequence, normalize=True)

    def sliding_window_analysis(self, sequence: str) -> List[Dict]:
        """Perform sliding window analysis on longer sequence"""
        return self.scorer.sliding_window_score(sequence, normalize=True)

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

                results.append(
                    {
                        "gene_id": region["gene_id"],
                        "upstream_sequence": sequence,
                        "best_score": score,
                        "best_position": 0,
                        "best_subsequence": sequence,
                        "has_promoter": score > self.threshold,
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
