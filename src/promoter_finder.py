"""
Manual promoter extraction based on WAWWWT pattern
Following Liu et al. 2011 methodology from assignment
"""

import logging
import re
from typing import Dict, List, Optional


class PromoterFinder:
    def __init__(self, filter_mode: str = "consecutive_w"):
        """Initialize promoter finder for WAWWWT pattern
        
        Args:
            filter_mode: "consecutive_w" for WWWWWW pattern (default)
                        "wawwwt" for strict WAWWWT pattern matching
        """
        self.filter_mode = filter_mode
        self.w_pattern = "[AT]"
        logging.info(f"Initialized PromoterFinder with filter_mode={filter_mode}")

    def has_consecutive_w_bases(self, sequence: str, min_length: int = 6) -> bool:
        """Check if sequence has at least 6 consecutive W (A or T) bases"""
        consecutive_w = re.search(f"[AT]{{{min_length},}}", sequence)
        return consecutive_w is not None
    
    def has_wawwwt_pattern(self, sequence: str) -> bool:
        """Check if sequence contains WAWWWT pattern: W-A-W-W-W-T"""
        wawwwt_pattern = re.search(r"[AT]A[AT][AT][AT]T", sequence)
        return wawwwt_pattern is not None

    def find_promoter_in_sequence(self, sequence: str) -> Optional[Dict]:
        """Find best promoter match in a sequence based on filter mode"""
        if self.filter_mode == "consecutive_w":
            if not self.has_consecutive_w_bases(sequence):
                return None
        elif self.filter_mode == "wawwwt":
            if not self.has_wawwwt_pattern(sequence):
                return None

        best_match = None
        best_score = 0

        for i in range(len(sequence) - 5):
            subseq = sequence[i : i + 6]

            if self.filter_mode == "wawwwt":
                if re.match(r"[AT]A[AT][AT][AT]T", subseq):
                    w_count = sum(1 for base in subseq if base in "AT")
                else:
                    continue
            else:
                w_count = sum(1 for base in subseq if base in "AT")

            if w_count > best_score:
                best_score = w_count
                best_match = {
                    "sequence": subseq,
                    "position": i,
                    "w_content": w_count,
                    "score": best_score,
                }

        return best_match

    def extract_promoters_manual(
        self, upstream_regions: List[Dict], sample_size: int = 100
    ) -> List[Dict]:
        """Manually extract promoters from sequences with best W content"""
        promoters = []
        candidates = []

        logging.info(
            f"Starting manual promoter extraction from up to {len(upstream_regions)} sequences (filter_mode={self.filter_mode})"
        )

        for region in upstream_regions:
            sequence = region["upstream_sequence"]

            if self.filter_mode == "consecutive_w":
                if not self.has_consecutive_w_bases(sequence):
                    continue
            elif self.filter_mode == "wawwwt":
                if not self.has_wawwwt_pattern(sequence):
                    continue

            promoter_match = self.find_promoter_in_sequence(sequence)

            if promoter_match:
                candidates.append(
                    {
                        "gene_id": region["gene_id"],
                        "promoter_sequence": promoter_match["sequence"],
                        "position_in_upstream": promoter_match["position"],
                        "w_content": promoter_match["w_content"],
                        "score": promoter_match["score"],
                        "original_sequence": sequence,
                    }
                )

        candidates.sort(key=lambda x: x["score"], reverse=True)
        promoters = candidates[:sample_size]

        filter_desc = "6 consecutive Ws" if self.filter_mode == "consecutive_w" else "WAWWWT pattern"
        logging.info(f"Found {len(candidates)} candidates with {filter_desc}")
        logging.info(f"Selected top {len(promoters)} promoters for training")
        return promoters
