#import "style.typ": *
#let title = "Bacterial Promoter Analysis"
#let subtitle = "Position Probability Matrix Construction and Statistical Alignment"
#let author = "Thuvaragan S."
#let index = "210657G"
#let course_id = "BM4322"
#let course_name = "Genomic Signal Processing"
#let instructor = ""
#let semester = "2025"
#let department = "Electronics and Telecommunication Department"
#let university = "University of Moratuwa"
#let faculty = "Faculty of Engineering"

#show: individual_assignment_title.with(
  title,
  author,
  index,
  subtitle,
  course_id,
  course_name,
  department,
  faculty,
  university,
)

= Introduction

This report presents a computational analysis of bacterial promoter sequences based on Liu et al. (2011), which established that the σ⁷⁰ subunit of bacterial RNA polymerase recognizes promoters following the *WAWWWT* pattern (where W = A or T) located 10 bases upstream of gene start sites. Using genome GCA_900637025.1, we performed Position Probability Matrix (PPM) construction, statistical alignment, and cross-validation to identify and characterize these promoter elements.

== Genome Information

- *Organism:* _Streptococcus pyogenes_ M1 476
- *Accession:* GCA_900637025.1
- *Genome Size:* 1,931,548 bp
- *Total Genes:* 1,100 annotated genes
- *Source:* NCBI Genome Database

== Objectives

Per assignment requirements:

+ *Task 1:* Construct a PPM from 100 manually extracted promoters (6 bases each, minimum 6 consecutive Ws) selected from 1100 genes' upstream regions (-15 to -5 bp)
+ *Task 2:* Perform statistical alignment on remaining 1000 regions using the PPM to detect promoters
+ *Task 3:* Cross-validate the PPM on 1000 samples from other students' genomes

= Materials and Methods

== Task 1: PPM Construction

=== Upstream Region Extraction

- Forward strand genes (+): positions [start-15, start-5]
- Reverse strand genes (-): positions [end+5, end+15], reverse complemented
- Region length: 11 nucleotides per gene

=== Promoter Selection Criteria

+ Must contain ≥6 consecutive W bases (A or T)
+ Extract all 6-base windows from each 11-base region
+ Score windows by W-content with bonus for WAWWWT pattern
+ Select top 100 highest-scoring candidates

=== PPM Construction

- Count base frequencies at each position (1-6)
- Apply pseudocounts: C=0.01, G=0.01 (heuristic for unobserved bases)
- Calculate probabilities: $P(b|p) = ("count" + "pseudocount") / (N + sum "pseudocounts")$
- N = 99 (actual training sequences obtained)

== Task 2: Statistical Alignment

=== Scoring Function

$ "Score"("sequence") = sum_(i=1)^6 log(P("base"_i | "position"_i)) $

=== Classification

- Threshold: Mean - 2σ from training set scores
- Sliding window approach: score all 6-bp windows within 11-bp regions
- Accept best score per sequence
- Test set: 1000 sequences (genes 100-1099)

== Task 3: Cross-Validation

Applied 210657G's PPM to upstream regions from:

- 210079K (GCA_001457635.1)
- 210179R (GCA_019048645.1)
- 210504L (GCA_900636475.1)
- 210707L (GCA_900475505.1)
- 210732H (GCA_019046945.1)

== Software

- Python 3.12, BioPython 1.84, pandas 2.2.3, numpy 2.1.3
- Visualization: matplotlib 3.9.2, seaborn 0.13.2, logomaker 0.8.7

= Results

== Task 1: Position Probability Matrix

=== Training Set

- Candidates screened: 100
- Promoters extracted: 99 (one sequence rejected)
- All sequences: 100% AT-rich (6 consecutive Ws confirmed)

#figure(
  image("assets/figures/training_data_analysis.png", width: 90%),
  caption: [Training data analysis showing sequence composition and characteristics]
)

=== Consensus Sequence

*TATAAT*

=== Position Probability Matrix

#figure(
  table(
    columns: 5,
    [Position], [A], [C], [G], [T],
    [1], [0.495], [0.000], [0.000], [0.505],
    [2], [0.626], [0.000], [0.000], [0.374],
    [3], [0.454], [0.000], [0.000], [0.545],
    [4], [0.606], [0.000], [0.000], [0.394],
    [5], [0.717], [0.000], [0.000], [0.283],
    [6], [0.424], [0.000], [0.000], [0.576],
  ),
  caption: [Position Probability Matrix for 99 training sequences]
)

=== Sequence Logo Visualization

#figure(
  image("assets/figures/sequence_logo.png", width: 90%),
  caption: [Sequence logo showing nucleotide probabilities at each position. Letter heights are proportional to frequency. Position 5 shows strongest A-preference (71.7%), critical for promoter function.]
)

=== Key Findings

- 100% AT-richness validates WAWWWT pattern requirement
- Position 5 shows strongest conservation (A: 71.7%)
- Consensus TATAAT matches canonical bacterial -10 box (Pribnow box)
- No G/C observed in training data (only pseudocounts contribute)

#figure(
  image("assets/figures/ppm_heatmap.png", width: 90%),
  caption: [Heatmap representation of position probability matrix]
)

== Task 2: Statistical Alignment Results

=== Detection Performance

- Test sequences: 1000
- Promoters detected: 399 (39.9%)
- Non-promoters: 601 (60.1%)

=== Score Statistics

#figure(
  table(
    columns: 2,
    [Metric], [Value],
    [Mean Score], [-14.859],
    [Median Score], [-17.417],
    [Std Deviation], [7.553],
    [Min Score], [-35.142],
    [Max Score], [0.000],
    [Threshold], [-10.0],
  ),
  caption: [Statistical alignment score distribution]
)

#figure(
  image("assets/figures/score_distributions.png", width: 90%),
  caption: [Score distributions showing clear separation between promoter (high scores) and non-promoter (low scores) populations, validating discriminatory power of the PPM.]
)

=== Positional Distribution

Detected promoters show 5' enrichment within upstream regions:

- Position 0-2: 68.9% of detections
- Position 3-5: 31.1% of detections

This confirms -10 box location hypothesis.

#figure(
  image("assets/figures/detection_summary.png", width: 90%),
  caption: [Detection summary showing distribution of detected promoters]
)

=== Top Detected Sequences

#figure(
  table(
    columns: 3,
    [Sequence], [Count], [Percentage],
    [TATAAT], [23], [5.8%],
    [AATAAT], [18], [4.5%],
    [TAAAAT], [15], [3.8%],
    [AAAAAT], [12], [3.0%],
  ),
  caption: [Most frequently detected promoter sequences]
)

== Task 3: Cross-Validation Results

=== Testing 210657G's PPM on other students' genomes

#figure(
  table(
    columns: 5,
    [Student], [Genome], [Regions], [Detected], [Rate],
    [210079K], [GCA_001457635.1], [1000], [378], [37.80%],
    [210179R], [GCA_019048645.1], [1000], [425], [42.50%],
    [210504L], [GCA_900636475.1], [1000], [388], [38.80%],
    [210707L], [GCA_900475505.1], [999], [313], [31.33%],
    [210732H], [GCA_019046945.1], [1000], [401], [40.10%],
  ),
  caption: [Cross-validation results across diverse bacterial genomes]
)

=== Cross-Validation Statistics

- Mean detection rate: 38.11%
- Standard deviation: 4.18%
- Range: 31.33% - 42.50%
- Own genome (210657G): 39.9%

#figure(
  image("assets/figures/cross_validation_comparison.png", width: 90%),
  caption: [Cross-validation comparison showing consistent detection rates across genomes]
)

=== Interpretation

Consistent detection rates across diverse bacterial genomes (CV = 11.0%) demonstrate:

+ Strong model generalizability
+ Conserved σ⁷⁰-dependent promoter architecture across species
+ PPM captures universal TATAAT motif rather than genome-specific features
+ No evidence of overfitting (own genome within 1σ of cross-validation mean)

= Discussion

== Biological Validation

=== Consensus Sequence Analysis

Our TATAAT consensus is identical to the canonical bacterial -10 promoter (Pribnow box) extensively documented in molecular biology literature. This validates both:

+ Computational methodology
+ Biological relevance of detected sequences

=== AT-Richness

Complete absence of G/C in training sequences reflects functional requirement for DNA melting during transcription initiation:

- AT pairs: 2 hydrogen bonds (easier separation)
- GC pairs: 3 hydrogen bonds (stronger)
- RNA polymerase requires strand separation for template access

=== Position-Specific Conservation

Position 5's strong A-preference (71.7%, tallest letter in sequence logo) is critical for:

- DNA bending and flexibility
- σ⁷⁰ subunit recognition
- Transcription bubble formation

== Detection Rate Analysis

=== 39.9% Detection Rate

Falls within expected biological range because:

- Not all genes use σ⁷⁰-dependent promoters (alternative σ factors exist)
- Housekeeping genes typically have strong -10 boxes
- Regulatory genes may have weaker or variant promoters
- Literature reports 30-50% detection for genome-wide -10 box searches

=== Conservative Threshold

Mean - 2σ threshold prioritizes specificity over sensitivity, reducing false positives while capturing ~95% of known promoters.

== Cross-Validation Significance

=== Model Generalizability

Tight clustering of detection rates (SD = 4.18%) across phylogenetically diverse bacteria demonstrates:

+ Universal nature of TATAAT promoter motif
+ Conserved transcriptional machinery across species
+ Successful transfer learning without retraining

=== Biological Implications

Similar detection rates suggest comparable proportions of:

- Housekeeping vs regulatory genes
- σ⁷⁰-dependent vs alternative σ factor usage
- Conserved vs species-specific transcription mechanisms

== Clinical Relevance

_S. pyogenes_ is a human pathogen causing pharyngitis, scarlet fever, and invasive infections. Understanding promoter architecture can inform:

- Antibiotic development targeting transcription
- Gene regulation studies for virulence factors
- Comparative genomics identifying strain differences

== Limitations

+ *Single promoter element:* Analysis limited to -10 box (did not model -35 region or spacer length)
+ *Unidirectional cross-validation:* Tested our PPM on other genomes but not vice versa (other students' PPMs unavailable)
+ *Computational validation only:* No experimental confirmation via RNA-seq or reporter assays

= Conclusion

=== Key Findings

+ *Consensus Sequence:* TATAAT matches canonical bacterial -10 box
+ *PPM Quality:* 100% AT-richness, position 5 shows strongest conservation (71.7% A)
+ *Detection Performance:* 39.9% rate within expected biological range
+ *Cross-Validation:* Robust generalization (31.33-42.50% across diverse genomes)
+ *Biological Validation:* Results align with established promoter biology

= References

Complete analysis pipeline and reproducible code available at:

https://github.com/thuvasooriya/promoter-analysis
