"""
Position Probability Matrix builder
Based on BM4321 statistical alignment lectures
"""

import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from scoring import SequenceScorer


class PPMBuilder:
    def __init__(self, pseudocount: float = 0.01):
        """Initialize PPM builder with pseudocount for all bases"""
        self.pseudocount = pseudocount
        self.bases = ["A", "C", "G", "T"]
        self.scorer: Optional[SequenceScorer] = None
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
            for base_idx in range(4):
                freq = frequency_matrix[base_idx, pos]
                ppm_matrix[base_idx, pos] = (freq + self.pseudocount) / (
                    num_sequences + 4 * self.pseudocount
                )

        # Create DataFrame for easy handling
        ppm_df = pd.DataFrame(
            ppm_matrix.T,
            columns=self.bases,
            index=[f"Position_{i + 1}" for i in range(seq_length)],
        )

        self.scorer = SequenceScorer(ppm_df)
        logging.info(f"Built PPM from {num_sequences} sequences of length {seq_length}")
        return ppm_df

    def get_consensus_sequence(self, ppm_df: pd.DataFrame) -> str:
        """Get consensus sequence from PPM"""
        if self.scorer is None:
            self.scorer = SequenceScorer(ppm_df)
        return self.scorer.get_consensus_sequence()

    def calculate_consensus_score(self, ppm_df: pd.DataFrame) -> float:
        """Calculate consensus sequence score"""
        if self.scorer is None:
            self.scorer = SequenceScorer(ppm_df)
        return self.scorer.calculate_consensus_score()

    def score_sequence(self, sequence: str, ppm_df: pd.DataFrame) -> float:
        """Score a sequence using PPM"""
        if self.scorer is None:
            self.scorer = SequenceScorer(ppm_df)
        return self.scorer.score_sequence(sequence, normalize=True)
