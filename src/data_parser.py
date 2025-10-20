"""
Updated data parsing utilities for genome analysis
Adapted for your specific data structure
"""

import logging
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO


class GenomeParser:
    def __init__(self, genome_id="GCA_900637025.1", student_id="210657G"):
        """Initialize parser with genome ID and student ID"""
        self.genome_id = genome_id
        self.student_id = student_id
        self.data_dir = Path("data") / f"{student_id}_{genome_id}" / genome_id

        self.fasta_file = self._find_fasta_file()
        self.gff_file = self.data_dir / "genomic.gff"

        self.genome_record = None
        self.genes = []

        logging.info(f"Initialized parser for student {student_id}, genome {genome_id}")
        logging.info(f"FASTA file: {self.fasta_file}")
        logging.info(f"GFF file: {self.gff_file}")

    def _find_fasta_file(self):
        """Find the FASTA file with correct naming pattern"""
        patterns = [
            f"{self.genome_id}_46338_H01_genomic.fna",
            f"{self.genome_id}_NCTC7465_genomic.fna",
            f"{self.genome_id}_ASM1904864v1_genomic.fna",
            f"{self.genome_id}_42197_F01_genomic.fna",
            f"{self.genome_id}_42925_G01_genomic.fna",
            f"{self.genome_id}_ASM1904694v1_genomic.fna",
        ]

        for pattern in patterns:
            candidate = self.data_dir / pattern
            if candidate.exists():
                return candidate

        fna_files = list(self.data_dir.glob("*.fna"))
        if fna_files:
            return fna_files[0]

        raise FileNotFoundError(f"No FASTA file found in {self.data_dir}")

    def load_genome(self):
        """Load genome sequence from FASTA file"""
        try:
            if not self.fasta_file.exists():
                raise FileNotFoundError(f"FASTA file not found: {self.fasta_file}")

            # Handle multi-contig genomes
            records = list(SeqIO.parse(self.fasta_file, "fasta"))

            if len(records) == 1:
                self.genome_record = records[0]
                logging.info(
                    f"Loaded single contig genome of length {len(self.genome_record.seq)}"
                )
            else:
                # Concatenate contigs for simplicity (or handle separately if needed)
                self.genome_records = records
                self.genome_record = records[0]  # Use first contig as primary
                logging.info(
                    f"Loaded {len(records)} contigs, using first contig (length: {len(self.genome_record.seq)})"
                )

            return self.genome_record

        except Exception as e:
            logging.error(f"Error loading genome: {e}")
            raise

    def parse_genes_from_gff(self, max_genes=1100):
        """Parse gene annotations from GFF file"""
        genes = []

        try:
            if not self.gff_file.exists():
                raise FileNotFoundError(f"GFF file not found: {self.gff_file}")

            with open(self.gff_file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith("#") or line.strip() == "":
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 9:
                        continue

                    # Look for gene features
                    feature_type = parts[2].lower()
                    if feature_type == "gene":
                        try:
                            gene_info = {
                                "chromosome": parts[0],
                                "source": parts[1],
                                "feature": parts[2],
                                "start": int(parts[3])
                                - 1,  # Convert to 0-based indexing
                                "end": int(parts[4]),
                                "score": parts[5] if parts[5] != "." else None,
                                "strand": parts[6],
                                "frame": parts[7] if parts[7] != "." else None,
                                "attributes": parts[8],
                                "line_number": line_num,
                            }

                            # Extract gene ID from attributes
                            attributes = parts[8]
                            gene_id_match = re.search(r"ID=([^;]+)", attributes)
                            if gene_id_match:
                                gene_info["gene_id"] = gene_id_match.group(1)
                            else:
                                gene_info["gene_id"] = f"gene_{len(genes) + 1}"

                            genes.append(gene_info)

                            if len(genes) >= max_genes:
                                break

                        except ValueError as ve:
                            logging.warning(f"Error parsing line {line_num}: {ve}")
                            continue

            self.genes = genes
            logging.info(f"Parsed {len(self.genes)} genes from GFF file")
            return self.genes

        except Exception as e:
            logging.error(f"Error parsing GFF file: {e}")
            raise

    def parse_genes(self, max_genes=1100):
        """Parse genes from GFF file"""
        return self.parse_genes_from_gff(max_genes)

    def extract_upstream_regions(self, upstream_start=-15, upstream_end=-5):
        """Extract upstream regions for all genes"""
        upstream_regions = []

        if not self.genome_record:
            self.load_genome()

        if not self.genes:
            self.parse_genes()

        if not self.genome_record:
            raise RuntimeError("Failed to load genome record")

        genome_length = len(self.genome_record.seq)

        for i, gene in enumerate(self.genes):
            try:
                if gene["strand"] == "+":
                    # Forward strand: upstream is before start
                    region_start = gene["start"] + upstream_start  # -15 to -5
                    region_end = gene["start"] + upstream_end + 1  # Include -5 position
                else:
                    # Reverse strand: upstream is after end
                    region_start = gene["end"] - upstream_end - 1  # Include +5 position
                    region_end = gene["end"] - upstream_start  # +15 position

                # Ensure coordinates are within genome bounds
                region_start = max(0, region_start)
                region_end = min(genome_length, region_end)

                if region_start < region_end:
                    upstream_seq = str(self.genome_record.seq[region_start:region_end])

                    if gene["strand"] == "-":
                        # Reverse complement for reverse strand
                        from Bio.Seq import Seq

                        upstream_seq = str(Seq(upstream_seq).reverse_complement())

                    upstream_regions.append(
                        {
                            "gene_index": i,
                            "gene_id": gene["gene_id"],
                            "upstream_sequence": upstream_seq.upper(),
                            "genomic_start": region_start,
                            "genomic_end": region_end,
                            "strand": gene["strand"],
                            "gene_start": gene["start"],
                            "gene_end": gene["end"],
                            "sequence_length": len(upstream_seq),
                        }
                    )
                else:
                    logging.warning(
                        f"Invalid coordinates for gene {gene['gene_id']}: start={region_start}, end={region_end}"
                    )

            except Exception as e:
                logging.warning(
                    f"Error processing gene {gene.get('gene_id', 'unknown')}: {e}"
                )
                continue

        logging.info(f"Extracted {len(upstream_regions)} upstream regions")
        logging.info(
            f"Average upstream sequence length: {np.mean([r['sequence_length'] for r in upstream_regions]):.1f}"
        )

        return upstream_regions
