"""
Centralized scoring utilities for promoter analysis
"""

import logging
from typing import List

import numpy as np
import pandas as pd


class SequenceScorer:
    def __init__(self, ppm_df: pd.DataFrame, filter_mode: str = "consecutive_w"):
        """Initialize scorer with Position Probability Matrix
        
        Args:
            ppm_df: Position Probability Matrix DataFrame
            filter_mode: "consecutive_w" or "wawwwt" (for logging/tracking)
        """
        self.ppm_df = ppm_df
        self.ppm_length = len(ppm_df)
        self.bases = ["A", "C", "G", "T"]
        self.filter_mode = filter_mode
        self.consensus_score = self.calculate_consensus_score()
        logging.info(f"Initialized SequenceScorer with PPM of length {self.ppm_length}, filter_mode={filter_mode}")

    def calculate_consensus_score(self) -> float:
        """Calculate consensus sequence benchmark score"""
        consensus_probs = []
        for _, row in self.ppm_df.iterrows():
            max_prob = row.max()
            consensus_probs.append(np.log(max_prob))
        return sum(consensus_probs)

    def get_consensus_sequence(self) -> str:
        """Get consensus sequence from PPM"""
        consensus = ""
        for _, row in self.ppm_df.iterrows():
            max_base = str(row.idxmax())
            consensus += max_base
        return consensus

    def score_sequence(self, sequence: str, normalize: bool = True) -> float:
        """
        Score a sequence using PPM
        
        Args:
            sequence: DNA sequence to score
            normalize: If True, subtract consensus score (default)
        
        Returns:
            Log-probability score (normalized if normalize=True)
        """
        if len(sequence) != self.ppm_length:
            raise ValueError(
                f"Sequence length {len(sequence)} != PPM length {self.ppm_length}"
            )

        log_score = 0.0
        for pos, base in enumerate(sequence):
            if base in self.bases:
                prob = self.ppm_df.iloc[pos][base]
                if prob > 0:
                    log_score += np.log(prob)
                else:
                    log_score += np.log(1e-10)
            else:
                logging.warning(f"Unknown base {base} at position {pos}")
                log_score += np.log(1e-10)

        if normalize:
            return log_score - self.consensus_score
        else:
            return log_score

    def score_sequences(self, sequences: List[str], normalize: bool = True) -> List[float]:
        """Score multiple sequences"""
        return [self.score_sequence(seq, normalize) for seq in sequences]

    def sliding_window_score(self, sequence: str, normalize: bool = True) -> List[dict]:
        """
        Perform sliding window analysis on longer sequence
        
        Returns:
            List of dicts with keys: position, sequence, score
        """
        results = []

        if len(sequence) < self.ppm_length:
            return results

        for i in range(len(sequence) - self.ppm_length + 1):
            subseq = sequence[i : i + self.ppm_length]
            score = self.score_sequence(subseq, normalize)

            results.append({"position": i, "sequence": subseq, "score": score})

        return results
