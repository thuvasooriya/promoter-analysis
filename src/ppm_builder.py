"""
Position Probability Matrix builder
Based on BM4321 statistical alignment lectures
"""

import logging
from typing import List

import numpy as np
import pandas as pd


class PPMBuilder:
    def __init__(self, pseudocount: float = 0.01):
        """Initialize PPM builder with pseudocount for C and G"""
        self.pseudocount = pseudocount
        self.bases = ["A", "C", "G", "T"]
        logging.info(f"Initialized PPM builder with pseudocount {pseudocount}")

    def build_ppm(self, promoter_sequences: List[str]) -> pd.DataFrame:
        """Build Position Probability Matrix from promoter sequences"""
        if not promoter_sequences:
            raise ValueError("No promoter sequences provided")

        seq_length = len(promoter_sequences[0])
        num_sequences = len(promoter_sequences)

        # Initialize frequency matrix
        frequency_matrix = np.zeros((4, seq_length))

        # Count base frequencies at each position
        for seq in promoter_sequences:
            for pos, base in enumerate(seq):
                if base in self.bases:
                    base_idx = self.bases.index(base)
                    frequency_matrix[base_idx, pos] += 1

        # Convert to probabilities with pseudocounts
        ppm_matrix = np.zeros((4, seq_length))

        for pos in range(seq_length):
            for base_idx, base in enumerate(self.bases):
                freq = frequency_matrix[base_idx, pos]

                # Add pseudocount for C and G (heuristic values)
                if base in ["C", "G"]:
                    freq += self.pseudocount

                # Normalize by total + pseudocounts
                total = num_sequences + (
                    2 * self.pseudocount
                )  # Two bases get pseudocount
                ppm_matrix[base_idx, pos] = freq / total

        # Create DataFrame for easy handling
        ppm_df = pd.DataFrame(
            ppm_matrix.T,
            columns=self.bases,
            index=[f"Position_{i + 1}" for i in range(seq_length)],
        )

        logging.info(f"Built PPM from {num_sequences} sequences of length {seq_length}")
        return ppm_df

    def get_consensus_sequence(self, ppm_df: pd.DataFrame) -> str:
        """Get consensus sequence from PPM"""
        consensus = ""
        for _, row in ppm_df.iterrows():
            max_base = row.idxmax()
            consensus += max_base
        return consensus

    def calculate_consensus_score(self, ppm_df: pd.DataFrame) -> float:
        """Calculate consensus sequence score"""
        consensus = self.get_consensus_sequence(ppm_df)
        return self.score_sequence(consensus, ppm_df)

    def score_sequence(self, sequence: str, ppm_df: pd.DataFrame) -> float:
        """Score a sequence using the PPM (log probabilities)"""
        if len(sequence) != len(ppm_df):
            raise ValueError(
                f"Sequence length {len(sequence)} != PPM length {len(ppm_df)}"
            )

        log_score = 0.0
        for pos, base in enumerate(sequence):
            if base in self.bases:
                prob = ppm_df.iloc[pos][base]
                log_score += np.log(prob)
            else:
                logging.warning(f"Unknown base {base} at position {pos}")

        return log_score
