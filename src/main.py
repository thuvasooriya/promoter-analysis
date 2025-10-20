#!/usr/bin/env python3
"""
BM4322 Genomic Signal Processing Assignment
Student: 210657G
Main analysis entry point
"""

import argparse
import json
import logging
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from data_parser import GenomeParser
from ppm_builder import PPMBuilder
from promoter_finder import PromoterFinder
from statistical_alignment import StatisticalAligner


class BM4322Analysis:
    def __init__(self, genome_id="GCA_900637025.1", student_id="210657G"):
        """Initialize analysis for specific genome and student"""
        self.genome_id = genome_id
        self.student_id = student_id

        self.data_dir = Path("data")
        self.results_dir = Path("results")
        self.task1_dir = Path("results/task1_ppm")
        self.task2_dir = Path("results/task2_alignment")
        self.task3_dir = Path("results/task3_cross_validation")
        self.logs_dir = Path("logs")

        for dir_path in [
            self.results_dir,
            self.task1_dir,
            self.task2_dir,
            self.task3_dir,
            self.logs_dir,
        ]:
            dir_path.mkdir(exist_ok=True, parents=True)

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler(self.logs_dir / "main_analysis.log"),
                logging.StreamHandler(),
            ],
        )

        logging.info(
            f"Initialized analysis for student {student_id}, genome {genome_id}"
        )

    def verify_data_files(self):
        """Verify that required data files exist"""
        parser = GenomeParser(self.genome_id, self.student_id)
        required_files = [parser.fasta_file, parser.gff_file]

        missing_files = [f for f in required_files if not f.exists()]

        if missing_files:
            logging.error("Missing required files:")
            for f in missing_files:
                logging.error(f"  - {f}")
            raise FileNotFoundError(f"Missing {len(missing_files)} required files")

        logging.info("All required data files found")
        return True

    def run_task1(self):
        """Run Task 1: Extract promoters and build PPM"""
        logging.info("=== Task 1: PPM Construction ===")

        self.verify_data_files()

        parser = GenomeParser(self.genome_id, self.student_id)
        genome_record = parser.load_genome()
        genes = parser.parse_genes(max_genes=1100)
        upstream_regions = parser.extract_upstream_regions()

        upstream_df = pd.DataFrame(upstream_regions)
        upstream_df.to_csv(self.task1_dir / "upstream_regions.csv", index=False)

        genes_df = pd.DataFrame(genes)
        genes_df.to_csv(self.task1_dir / "genes_info.csv", index=False)

        logging.info(f"Genome length: {len(genome_record.seq):,} bp")
        logging.info(f"Total genes found: {len(genes)}")
        logging.info(f"Upstream regions extracted: {len(upstream_regions)}")

        promoter_finder = PromoterFinder()
        valid_regions = [r for r in upstream_regions if r["sequence_length"] == 11]
        logging.info(f"Valid 11-bp upstream regions: {len(valid_regions)}")

        if len(valid_regions) < 100:
            logging.warning(f"Only {len(valid_regions)} valid regions, using all")
            sample_size = len(valid_regions)
        else:
            sample_size = 100

        promoters = promoter_finder.extract_promoters_manual(valid_regions, sample_size)

        promoters_df = pd.DataFrame(promoters)
        promoters_df.to_csv(self.task1_dir / "manual_promoters.csv", index=False)

        logging.info(
            f"Extracted {len(promoters)} promoters from {sample_size} sequences"
        )

        if len(promoters) < 10:
            logging.error("Too few promoters found!")
            return None, None, None

        ppm_builder = PPMBuilder(pseudocount=0.01)
        promoter_sequences = [p["promoter_sequence"] for p in promoters]
        ppm_df = ppm_builder.build_ppm(promoter_sequences)

        ppm_df.to_csv(self.task1_dir / "position_probability_matrix.csv")

        consensus = ppm_builder.get_consensus_sequence(ppm_df)
        consensus_score = ppm_builder.calculate_consensus_score(ppm_df)

        logging.info(f"Consensus sequence: {consensus}")
        logging.info(f"Consensus score: {consensus_score:.3f}")

        task1_summary = {
            "student_id": self.student_id,
            "genome_id": self.genome_id,
            "genome_length": len(genome_record.seq),
            "total_genes": len(genes),
            "upstream_regions": len(upstream_regions),
            "valid_regions": len(valid_regions),
            "training_size": sample_size,
            "promoters_found": len(promoters),
            "consensus_sequence": consensus,
            "consensus_score": float(consensus_score),
        }

        with open(self.task1_dir / "task1_summary.json", "w") as f:
            json.dump(task1_summary, f, indent=2)

        logging.info("=== Task 1 Complete ===")
        return ppm_df, valid_regions, sample_size

    def run_task2(self, ppm_df=None, valid_regions=None, sample_size=None):
        """Run Task 2: Statistical alignment on remaining regions"""
        logging.info("=== Task 2: Statistical Alignment ===")

        if ppm_df is None:
            ppm_file = self.task1_dir / "position_probability_matrix.csv"
            if not ppm_file.exists():
                logging.error("PPM not found. Run Task 1 first.")
                return
            ppm_df = pd.read_csv(ppm_file, index_col=0)

        if valid_regions is None or sample_size is None:
            parser = GenomeParser(self.genome_id, self.student_id)
            parser.load_genome()
            parser.parse_genes(max_genes=1100)
            upstream_regions = parser.extract_upstream_regions()
            valid_regions = [r for r in upstream_regions if r["sequence_length"] == 11]
            sample_size = 100

        aligner = StatisticalAligner(ppm_df)
        test_regions = valid_regions[sample_size : sample_size + 1000]

        if len(test_regions) < 1000:
            logging.warning(f"Only {len(test_regions)} test sequences available")

        logging.info(f"Testing on {len(test_regions)} regions")
        alignment_results = aligner.analyze_upstream_regions(test_regions)

        results_df = pd.DataFrame(
            [
                {
                    "gene_id": r["gene_id"],
                    "upstream_sequence": r["upstream_sequence"],
                    "best_score": r["best_score"],
                    "best_position": r["best_position"],
                    "best_subsequence": r["best_subsequence"],
                    "has_promoter": r["has_promoter"],
                }
                for r in alignment_results
            ]
        )
        results_df.to_csv(
            self.task2_dir / "statistical_alignment_results.csv", index=False
        )

        promoter_count = sum(1 for r in alignment_results if r["has_promoter"])
        detection_rate = (
            promoter_count / len(alignment_results) if alignment_results else 0
        )

        task2_summary = {
            "test_regions": len(test_regions),
            "promoters_detected": promoter_count,
            "detection_rate": float(detection_rate),
            "threshold": float(aligner.threshold),
        }

        with open(self.task2_dir / "task2_summary.json", "w") as f:
            json.dump(task2_summary, f, indent=2)

        logging.info(
            f"Detected {promoter_count}/{len(test_regions)} promoters ({detection_rate:.1%})"
        )
        logging.info("=== Task 2 Complete ===")
        return alignment_results

    def run_task3(self):
        """Run Task 3: Cross-validation with other students' data"""
        logging.info("=== Task 3: Cross-Validation ===")

        students = {
            "210079K": "GCA_001457635.1",
            "210179R": "GCA_019048645.1",
            "210504L": "GCA_900636475.1",
            "210657G": "GCA_900637025.1",
            "210707L": "GCA_900475505.1",
            "210732H": "GCA_019046945.1",
        }

        ppm_file = self.task1_dir / "position_probability_matrix.csv"
        if not ppm_file.exists():
            ppm_file = self.results_dir / "position_probability_matrix.csv"
        if not ppm_file.exists():
            logging.error("PPM file not found. Run Task 1 first.")
            return

        our_ppm = pd.read_csv(ppm_file, index_col=0)

        cross_val_results = {}

        for student_id, genome_id in students.items():
            if student_id == self.student_id:
                continue

            logging.info(f"Testing PPM on {student_id}'s data ({genome_id})")

            try:
                test_results = self._test_ppm_on_student_data(
                    student_id, genome_id, our_ppm
                )
                cross_val_results[student_id] = test_results

            except Exception as e:
                logging.error(f"Failed to test on {student_id}: {e}")
                cross_val_results[student_id] = {"error": str(e)}

        cv_file = self.task3_dir / "cross_validation_results.json"
        with open(cv_file, "w") as f:
            json.dump(cross_val_results, f, indent=2)

        summary_lines = [
            "TASK 3: CROSS-VALIDATION RESULTS",
            "=" * 60,
            f"Testing {self.student_id}'s PPM on other students' genomes",
            "",
        ]

        for student_id, results in cross_val_results.items():
            if "error" in results:
                summary_lines.append(f"{student_id}: ERROR - {results['error']}")
            else:
                summary_lines.append(f"{student_id} ({results['genome_id']}):")
                summary_lines.append(
                    f"  Regions tested: {results['total_regions_tested']}"
                )
                summary_lines.append(
                    f"  Promoters detected: {results['promoters_detected']}"
                )
                summary_lines.append(
                    f"  Detection rate: {results['detection_rate']:.2%}"
                )
                summary_lines.append("")

        summary_text = "\n".join(summary_lines)

        with open(self.task3_dir / "cross_validation_summary.txt", "w") as f:
            f.write(summary_text)

        logging.info("=== Task 3 Complete ===")

    def _test_ppm_on_student_data(self, target_student_id, target_genome_id, ppm):
        """Test PPM on another student's genome data"""
        temp_analysis = BM4322Analysis(target_genome_id, target_student_id)

        parser = GenomeParser(target_genome_id, target_student_id)
        genome_record = parser.load_genome()
        genes = parser.parse_genes(max_genes=1100)
        upstream_regions = parser.extract_upstream_regions()

        valid_regions = [r for r in upstream_regions if r["sequence_length"] == 11]
        test_regions = valid_regions[100:1100]

        aligner = StatisticalAligner(ppm)
        alignment_results = aligner.analyze_upstream_regions(test_regions)

        promoter_count = sum(1 for r in alignment_results if r["has_promoter"])
        detection_rate = (
            promoter_count / len(alignment_results) if alignment_results else 0
        )

        return {
            "genome_id": target_genome_id,
            "total_regions_tested": len(alignment_results),
            "promoters_detected": promoter_count,
            "detection_rate": detection_rate,
            "score_statistics": {
                "mean": float(np.mean([r["best_score"] for r in alignment_results])),
                "std": float(np.std([r["best_score"] for r in alignment_results])),
                "max": float(max([r["best_score"] for r in alignment_results])),
                "min": float(min([r["best_score"] for r in alignment_results])),
            },
        }


def verify_data_structure(student_id="210657G", genome_id="GCA_900637025.1"):
    """Verify required data files exist"""
    data_dir = Path("data") / f"{student_id}_{genome_id}" / genome_id
    required_files = [
        data_dir / f"{genome_id}_46338_H01_genomic.fna",
        data_dir / "genomic.gff",
    ]

    missing_files = [f for f in required_files if not f.exists()]
    if missing_files:
        print("Missing required data files:")
        for f in missing_files:
            print(f"  - {f}")
        print("\nExpected structure:")
        print(f"  data/{student_id}_{genome_id}/{genome_id}/")
        print(f"    ├── {genome_id}_46338_H01_genomic.fna")
        print("    └── genomic.gff")
        return False
    return True


def main():
    print("BM4322 Genomic Signal Processing Assignment")
    print("Student ID: 210657G")
    print("Genome: GCA_900637025.1")
    print("=" * 50)

    student_id = "210657G"
    genome_id = "GCA_900637025.1"

    if not verify_data_structure(student_id, genome_id):
        return 1

    parser = argparse.ArgumentParser(description="BM4322 Assignment Runner")
    parser.add_argument(
        "--task",
        choices=["1", "2", "3", "all"],
        default="all",
        help="Task to run: 1 (PPM), 2 (alignment), 3 (cross-validation), all (complete)",
    )
    args = parser.parse_args()

    try:
        analysis = BM4322Analysis()

        if args.task in ["1", "all"]:
            analysis.run_task1()
            print("Task 1 (PPM construction) completed!")

        if args.task in ["2", "all"]:
            analysis.run_task2()
            print("Task 2 (statistical alignment) completed!")

        if args.task in ["3", "all"]:
            analysis.run_task3()
            print("Task 3 (cross-validation) completed!")

        print("\n" + "=" * 50)
        print("All requested tasks completed!")
        print("Results saved in:")
        if args.task in ["1", "all"]:
            print("  - results/task1_ppm/")
        if args.task in ["2", "all"]:
            print("  - results/task2_alignment/")
        if args.task in ["3", "all"]:
            print("  - results/task3_cross_validation/")

        return 0

    except Exception as e:
        print(f"Analysis failed with error: {e}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
