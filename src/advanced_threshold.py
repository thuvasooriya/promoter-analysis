"""
Advanced threshold selection methods for promoter detection
Uses training promoter scores to set detection threshold
Note: log-probability scores where HIGHER is BETTER (less negative)
"""

import logging
from typing import Dict, List, Tuple

import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture


class AdvancedThresholdCalculator:
    def __init__(self):
        """Initialize advanced threshold calculator"""
        self.methods = {
            "mean_minus_2std": self._mean_minus_2std,
            "percentile": self._percentile,
            "otsu": self._otsu,
            "gmm": self._gaussian_mixture,
            "iqr": self._iqr_method,
            "mad": self._mad_method,
        }
        logging.info("Initialized AdvancedThresholdCalculator")

    def _mean_minus_2std(self, scores: np.ndarray, **kwargs) -> float:
        """Classic mean - 2*std (lower tail for log-scores where lower = worse)"""
        mean = np.mean(scores)
        std = np.std(scores)
        return float(mean - 2 * std)

    def _percentile(self, scores: np.ndarray, percentile: float = 5, **kwargs) -> float:
        """Percentile-based threshold (default 5th percentile = lower tail)"""
        return float(np.percentile(scores, percentile))

    def _otsu(self, scores: np.ndarray, **kwargs) -> float:
        """Otsu's method - finds threshold that minimizes intra-class variance"""
        scores_sorted = np.sort(scores)
        n = len(scores_sorted)

        if n < 10:
            return self._mean_minus_2std(scores)

        best_threshold = scores_sorted[0]
        best_variance = float("inf")

        for i in range(1, n - 1):
            threshold = scores_sorted[i]

            left = scores_sorted[:i]
            right = scores_sorted[i:]

            if len(left) < 2 or len(right) < 2:
                continue

            w1 = len(left) / n
            w2 = len(right) / n

            var1 = np.var(left)
            var2 = np.var(right)

            within_class_variance = w1 * var1 + w2 * var2

            if within_class_variance < best_variance:
                best_variance = within_class_variance
                best_threshold = threshold

        return best_threshold

    def _gaussian_mixture(self, scores: np.ndarray, **kwargs) -> float:
        """Gaussian Mixture Model - assumes bimodal distribution, finds lower boundary"""
        if len(scores) < 20:
            return self._mean_minus_2std(scores)

        scores_reshaped = scores.reshape(-1, 1)

        try:
            gmm = GaussianMixture(n_components=2, random_state=42)
            gmm.fit(scores_reshaped)

            means = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_.flatten())

            low_idx = np.argmin(means)

            threshold = float(means[low_idx] + 2 * stds[low_idx])

            return threshold

        except Exception as e:
            logging.warning(f"GMM failed: {e}, falling back to mean-2std")
            return self._mean_minus_2std(scores)

    def _iqr_method(self, scores: np.ndarray, **kwargs) -> float:
        """Interquartile Range method (lower outlier detection)"""
        q1 = np.percentile(scores, 25)
        q3 = np.percentile(scores, 75)
        iqr = q3 - q1

        threshold = float(q1 - 1.5 * iqr)
        return threshold

    def _mad_method(self, scores: np.ndarray, **kwargs) -> float:
        """Median Absolute Deviation (very robust, lower tail)"""
        median = np.median(scores)
        mad = np.median(np.abs(scores - median))

        threshold = float(median - 2.5 * mad)
        return threshold

    def calculate_threshold(
        self, training_scores: List[float], method: str = "mean_minus_2std", **kwargs
    ) -> Tuple[float, Dict]:
        """
        Calculate threshold using specified method

        Args:
            training_scores: List of scores from training promoters
            method: One of the available methods
            **kwargs: Method-specific parameters

        Returns:
            Tuple of (threshold, statistics_dict)
        """
        if method not in self.methods:
            available = ", ".join(self.methods.keys())
            raise ValueError(f"Unknown method '{method}'. Available: {available}")

        scores = np.array(training_scores)

        threshold = self.methods[method](scores, **kwargs)

        stats_dict = {
            "method": method,
            "threshold": float(threshold),
            "training_mean": float(np.mean(scores)),
            "training_median": float(np.median(scores)),
            "training_std": float(np.std(scores)),
            "training_min": float(np.min(scores)),
            "training_max": float(np.max(scores)),
            "training_count": len(scores),
        }

        logging.info(
            f"Threshold calculated using {method}: {threshold:.3f} "
            f"(training mean: {stats_dict['training_mean']:.3f}, "
            f"std: {stats_dict['training_std']:.3f})"
        )

        return threshold, stats_dict

    def compare_methods(
        self, training_scores: List[float]
    ) -> Dict[str, Tuple[float, Dict]]:
        """
        Compare all threshold calculation methods

        Returns:
            Dictionary mapping method name to (threshold, stats_dict)
        """
        results = {}

        for method_name in self.methods.keys():
            try:
                threshold, stats = self.calculate_threshold(
                    training_scores, method=method_name
                )
                results[method_name] = (threshold, stats)
            except Exception as e:
                logging.error(f"Method {method_name} failed: {e}")
                results[method_name] = (None, {"error": str(e)})

        return results
