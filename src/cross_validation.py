"""
Cross-validation module for BM4322 assignment Task 3
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from main import BM4322Analysis
from statistical_alignment import StatisticalAligner


class CrossValidator:
    """Handles cross-validation between students' PPMs and test data"""

    def __init__(self):
        self.students = {
            "210079K": "GCA_001457635.1",
            "210179R": "GCA_019048645.1",
            "210504L": "GCA_900636475.1",
            "210657G": "GCA_900637025.1",
            "210707L": "GCA_900475505.1",
            "210732H": "GCA_019046945.1",
        }
        self.results_dir = Path("results")
        self.results_dir.mkdir(exist_ok=True)

    def run_complete_cross_validation(self):
        """Run all cross-validations and generate comprehensive report"""
        logging.info("Starting complete cross-validation analysis")

        all_results = {}

        for train_student, train_genome in self.students.items():
            logging.info(f"Loading PPM from {train_student}")

            # Load trained PPM
            ppm = self._load_student_ppm(train_student)
            if ppm is None:
                logging.warning(f"No PPM found for {train_student}, skipping")
                continue

            student_results = {}

            # Test on all other students
            for test_student, test_genome in self.students.items():
                if test_student == train_student:
                    continue

                logging.info(f"Testing {train_student}'s PPM on {test_student}'s data")
                result = self._validate_on_student(test_student, test_genome, ppm)
                student_results[test_student] = result

            all_results[train_student] = student_results

        # Save comprehensive results
        output_file = self.results_dir / "complete_cross_validation.json"
        with open(output_file, "w") as f:
            json.dump(all_results, f, indent=2)

        # Generate summary report
        self._generate_summary_report(all_results)

        logging.info("Cross-validation analysis completed")
        return all_results

    def _load_student_ppm(self, student_id):
        """Load PPM from a student's results"""
        ppm_file = self.results_dir / f"{student_id}_position_probability_matrix.csv"
        if ppm_file.exists():
            return pd.read_csv(ppm_file, index_col=0)
        return None

    def _validate_on_student(self, student_id, genome_id, ppm):
        """Test PPM on a student's genome data"""
        try:
            # Create analysis instance for this student
            analysis = BM4322Analysis(genome_id, student_id)

            # Get test regions
            from data_parser import GenomeParser

            parser = GenomeParser(genome_id, student_id)

            genome_record = parser.load_genome()
            genes = parser.parse_genes(max_genes=1100)
            upstream_regions = parser.extract_upstream_regions()

            # Use all valid regions for testing
            valid_regions = [r for r in upstream_regions if r["sequence_length"] == 11]

            # Apply PPM
            aligner = StatisticalAligner(ppm)
            results = aligner.analyze_upstream_regions(valid_regions)

            # Calculate metrics
            promoter_count = sum(1 for r in results if r["has_promoter"])
            detection_rate = promoter_count / len(results) if results else 0

            return {
                "success": True,
                "regions_tested": len(results),
                "promoters_detected": promoter_count,
                "detection_rate": detection_rate,
                "mean_score": float(np.mean([r["best_score"] for r in results])),
            }

        except Exception as e:
            return {"success": False, "error": str(e)}

    def _generate_summary_report(self, results):
        """Generate human-readable summary of cross-validation results"""
        summary_lines = [
            "=" * 70,
            "BM4322 CROSS-VALIDATION RESULTS (Task 3)",
            "=" * 70,
            "",
            "Format: Train_PPM -> Test_Data : Detection_Rate",
            "",
        ]

        for train_student, test_results in results.items():
            summary_lines.append(f"PPM from {train_student}:")

            for test_student, result in test_results.items():
                if result.get("success", False):
                    rate = result["detection_rate"]
                    summary_lines.append(f"  -> {test_student}: {rate:.3f}")
                else:
                    summary_lines.append(f"  -> {test_student}: FAILED")

            summary_lines.append("")

        # Calculate overall statistics
        successful_tests = []
        for train_results in results.values():
            for test_result in train_results.values():
                if test_result.get("success", False):
                    successful_tests.append(test_result["detection_rate"])

        if successful_tests:
            summary_lines.extend(
                [
                    "OVERALL STATISTICS:",
                    f"  Mean detection rate: {np.mean(successful_tests):.3f}",
                    f"  Std detection rate: {np.std(successful_tests):.3f}",
                    f"  Min detection rate: {min(successful_tests):.3f}",
                    f"  Max detection rate: {max(successful_tests):.3f}",
                    f"  Successful cross-validations: {len(successful_tests)}/{len(successful_tests) * 6}",
                ]
            )

        summary_text = "\n".join(summary_lines)

        with open(self.results_dir / "cross_validation_summary.txt", "w") as f:
            f.write(summary_text)

