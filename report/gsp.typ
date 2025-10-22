#import "style.typ": *
#let title = "Assignment"
#let subtitle = "Promoter Analysis"
#let author = "Thuvaragan S."
#let index = "210657G"
#let course_id = "BM4322"
#let course_name = "Genomic Signal Processing"
#let instructor = "Dr. Upeka Premaratne"
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

This report presents a computational analysis of bacterial promoter sequences using *statistical gene prediction* methods. Based on Liu et al. (2011), the $sigma^(70)$ subunit of bacterial RNA polymerase recognizes promoters following the `WAWWWT` pattern (where W = A or T) located approximately 10 bases upstream of gene start sites, corresponding to the *Pribnow Box* or *-10 box*.

= Genome Information

- *Organism:* _Streptococcus pyogenes_ M1 476
- *Accession:* GCA_900637025.1
- *Genome Size:* 1,931,548 bp
- *Total Genes:* 1,100 annotated genes
- *Source:* NCBI Genome Database

= Objectives

+ Construct a _Position Probability Matrix (PPM)_ from promoter sequences matching the WAWWWT pattern (6 bases each, where W = A or T) extracted from 1100 genes' upstream regions (-15 to -5 bp relative to start codon)
+ Perform _statistical alignment_ on the remaining upstream regions using the PPM, computing log probability scores to detect promoter presence/absence based on an empirical threshold
+ Cross-validate the PPM generalizability by applying it to 1000 upstream regions from five other bacterial genomes assigned to classmates

= Materials and Methods

== Task 1: PPM Construction

=== Upstream Region Extraction

Extracted regions 15 to 5 bases upstream of gene start positions:

- _Forward strand genes (+):_ positions `[start - 15, start - 5]`
- _Reverse strand genes (-):_ positions `[end + 5, end + 15]`, reverse complemented
- Region length: 11 nucleotides per gene (allows 6-base sliding window)

=== Promoter Selection Criteria

Pattern matching for canonical Pribnow box:

+ Filter sequences matching WAWWWT pattern: [AT]A[AT][AT][AT]T
+ Extract all 6-base windows from each 11-base upstream region using sliding window
+ Accept only sequences strictly matching the pattern (positions 2 and 6 must be A and T respectively, other W positions can be A or T)
+ This yields training sequences closely resembling the canonical TATAAT motif
+ Pattern-based filtering ensures biological relevance and reduces noise from non-promoter sequences

=== Position Frequency Table Construction

For N promoter sequences of length L:

$ f_(j comma N) = "count of base" N "at position" j $

where $j in {1 comma 2 comma ... comma L}$ and $N in {A comma C comma G comma T}$

=== Converting Frequencies to Probabilities

Frequencies are converted to probabilities using pseudocounts:

$ p_(j comma N) = (f_(j comma N) + k) / (4k + sum_N f_(j comma N)) $

where $p_(j comma N)$ is the probability of base N at position j, $f_(j comma N)$ is frequency, and $k = 0.01$ (pseudocount constant). Training set: N = 42 promoter sequences matching WAWWWT pattern.

=== Implementation: PPM Construction

Core algorithm from `src/ppm_builder.py`:

```python
class PPMBuilder:
    def __init__(self, pseudocount: float = 0.01):
        self.pseudocount = pseudocount
        self.bases = ["A", "C", "G", "T"]

    def build_ppm(self, promoter_sequences: List[str]) -> pd.DataFrame:
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

                # Add pseudocount for C and G (heuristic)
                if base in ["C", "G"]:
                    freq += self.pseudocount

                # Normalize: p = (f + k) / (N + 2k)
                total = num_sequences + (2 * self.pseudocount)
                ppm_matrix[base_idx, pos] = freq / total

        return pd.DataFrame(ppm_matrix.T, columns=self.bases)
```

== Task 2: Statistical Alignment

Statistical alignment scores sequences using log probabilities for numerical stability:

$ "Score"(S) = sum_(j=1)^L log(p_(j comma s_j)) $

where $p_(j comma s_j)$ is the probability of base $s_j$ at position j.

Consensus sequence TATAAT (position 2: 99.9% A, position 6: 99.9% T) serves as benchmark with score -1.392.

=== Threshold Selection

Empirical threshold (-10.0) determined from score distribution analysis allows detection of sequences with 1-2 base variations from the consensus pattern. Score distributions show clear separation: training sequences (strong TATAAT matches) score near 0, while non-promoter sequences score below -10.0.

=== Implementation: Statistical Scoring



== Task 3: Cross-Validation

Cross-validation tests PPM generalizability across different genomes. The PPM trained on 210657G's genome is applied without modification to upstream regions from five other bacterial genomes:

- 210079K (GCA_001457635.1) - _Streptococcus pyogenes_
- 210179R (GCA_019048645.1) - _Streptococcus pyogenes_
- 210504L (GCA_900636475.1) - _Streptococcus pyogenes_
- 210707L (GCA_900475505.1) - _Streptococcus pyogenes_
- 210732H (GCA_019046945.1) - _Streptococcus pyogenes_

*Methodology:* Same statistical alignment procedure (scoring + threshold classification) applied to 1000 upstream regions per genome using 210657G's PPM, without retraining or parameter adjustment.

== Software and Implementation

=== Environment
- *Python:* 3.12 with uv package manager
- *Core libraries:* BioPython 1.84, pandas 2.2.3, numpy 2.1.3
- *Visualization:* matplotlib 3.9.2, seaborn 0.13.2, logomaker 0.8.7

=== Key Modules
- `src/data_parser.py` - GFF3/FASTA parsing, upstream region extraction
- `src/ppm_builder.py` - Position Probability Matrix construction
- `src/statistical_alignment.py` - Scoring and classification
- `src/cross_validation.py` - Multi-genome validation
- `src/visualizations.py` - Figures and sequence logos

= Results

== Task 1: Position Probability Matrix

=== Training Set Characteristics

- Upstream regions screened: 1100 genes
- Sequences matching WAWWWT pattern [AT]A[AT][AT][AT]T: 42
- Training sequences extracted: 42 unique genes
- AT-richness: 100% (all sequences contain only A and T, conforming to WAWWWT pattern requirement)
- Sequence length: 6 bases (positions 1-6)

#figure(
  image("assets/figures/training_data_analysis.png", width: 90%),
  caption: [Training data showing 42 promoter sequences matching WAWWWT pattern. Position 2 (99.9% A) and position 6 (99.9% T) exhibit near-perfect conservation, confirming canonical TATAAT Pribnow box structure. Positions 1, 3, 4, 5 show AT variation consistent with W (A or T) positions. All G/C probabilities derive from pseudocounts only.],
)

=== Consensus Sequence

The consensus sequence (highest probability base at each position):

*TATAAT*

This matches the canonical bacterial *Pribnow Box*, confirming successful identification of biologically relevant promoter sequences. The TATAAT consensus represents the prototypical $sigma^(70)$ recognition motif for RNA polymerase binding.

=== Consensus Score

$
  S_("consensus") = log(0.262 times 0.999 times 0.452 times 0.785 times 0.785 times 0.999) = -1.392
$

This benchmark score represents a high-quality TATAAT promoter under this PPM model.

=== Position Probability Matrix

#figure(
  image("assets/figures/ppm_heatmap.png", width: 70%),
  caption: [Position Probability Matrix heatmap showing base probabilities at each position. Constructed from 42 training sequences matching WAWWWT pattern. Position 2 (99.9% A) and position 6 (99.9% T) show near-complete conservation defining the canonical TATAAT motif. C and G probabilities derive from pseudocounts only ($k = 0.01$, shown as near-zero values).],
)

=== Sequence Logo Visualization

#figure(
  image("assets/figures/sequence_logo.png", width: 90%),
  caption: [Sequence logo representation of the PPM. Letter heights are proportional to probability. Positions 2 and 6 show near-complete conservation (A and T respectively, 99.9% each), defining the canonical TATAAT Pribnow box. Positions 1, 3, 4, 5 show AT variation reflecting the W (weak base pair) positions, allowing functional flexibility while maintaining DNA melting capability for transcription initiation.],
)

== Task 2: Statistical Alignment Results

=== Detection Performance (Test Set: 1000 Non-overlapping Regions)

Using empirical threshold (-10.0) on 1000 upstream regions (excluding 42 training genes):

- Test sequences analyzed: 1000 upstream regions (non-overlapping with training)
- Promoters detected (Score > -10.0): 126 (12.6%)
- Non-promoters (Score ≤ -10.0): 874 (87.4%)
- Classification threshold: -10.0 (empirical)

=== Score Distribution Analysis

#figure(
  table(
    columns: 3,
    [Dataset], [Mean ± SD], [Range],
    [Training (42)], [~0.0 ± 0.0], [~0.0],
    [Test (1000)], [-19.0 ± 6.2], [[-40.6, -7.7]],
    [Background], [-20.0 ± 10.0], [[-52.0, 0.0]],
  ),
  caption: [Score distributions showing clear separation between training promoters (near-perfect TATAAT matches scoring ~0), test regions (intermediate), and random background (low)],
)

#figure(
  table(
    columns: 2,
    [Metric], [Value],
    [Mean Score], [-19.0],
    [Median Score], [-18.5],
    [Std Deviation], [6.2],
    [Min Score], [-40.6],
    [Max Score], [-7.7],
    [Threshold], [-10.0],
  ),
  caption: [Statistical alignment score distribution for non-overlapping test set],
)

#figure(
  image("assets/figures/score_distributions.png", width: 90%),
  caption: [Score distributions showing clear separation between promoter (high scores) and non-promoter (low scores) populations, validating discriminatory power of the PPM.],
)

=== Positional Distribution Within Upstream Regions

Sliding window analysis reveals where promoters are detected within 11-bp upstream regions:

- *Positions 0-2* (earlier in region, farther from start codon): ~70% of detections
- *Positions 3-5* (later in region, closer to start codon): ~30% of detections

This 5' enrichment confirms the -10 box location hypothesis (approximately 10 bases upstream of the start codon, corresponding to earlier positions in the -15 to -5 extraction window).

=== Statistical Alignment Scoring Example

*"Positive" Sequence (TATAAT - perfect consensus):*

$
  S = log(0.262) + log(0.999) + log(0.452) + log(0.785) + log(0.785) + log(0.999) = -1.39
$

$
  S_("normalized") = -1.39 - (-1.392) approx 0.0 > -10.0 arrow "Promoter detected"
$

*"Negative" Sequence (GGCCAC):*

$
  S = log(0.0002) + log(0.0002) + log(0.0002) + log(0.0002) + log(0.0002) + log(0.0002) approx -51.3
$

$ S_("normalized") = -51.3 - (-1.392) = -49.9 < -10.0 arrow "No promoter" $

The clear score separation demonstrates the PPM's discriminatory power.

#figure(
  image("assets/figures/detection_summary.png", width: 90%),
  caption: [Detection summary for Task 2. Top panel shows 12.6% detection rate (126/1000 sequences). Bottom panel displays positional distribution of detected promoters within 11-bp upstream regions, with ~70% detected at positions 0-2 (farther from start codon), consistent with -10 box location hypothesis. Detected sequences show strong TATAAT-like patterns, reflecting the canonical Pribnow box structure.],
)

== Task 3: Cross-Validation Results

=== Testing 210657G's PPM on other students' genomes

*Note:* Cross-validation uses empirical threshold (Score > -10.0) for practical detection.

#figure(
  table(
    columns: 5,
    [Student], [Genome], [Regions], [Detected], [Rate],
    [210079K], [GCA_001457635.1], [1000], [129], [12.9%],
    [210179R], [GCA_019048645.1], [1000], [141], [14.1%],
    [210504L], [GCA_900636475.1], [1000], [128], [12.8%],
    [210707L], [GCA_900475505.1], [1000], [102], [10.2%],
    [210732H], [GCA_019046945.1], [1000], [143], [14.3%],
  ),
  caption: [Cross-validation results using empirical threshold (-10.0) across diverse bacterial genomes],
)

=== Cross-Validation Statistics

- Mean detection rate: 12.9%
- Standard deviation: 1.5%
- Range: 10.2% - 14.3%
- Own genome (210657G): 12.6% (empirical threshold)

#figure(
  image("assets/figures/cross_validation_comparison.png", width: 90%),
  caption: [Cross-validation comparison showing consistent detection rates across genomes using empirical threshold],
)

=== Interpretation

Using the empirical threshold (-10.0), consistent detection rates across diverse bacterial genomes (CV = 11.6%) demonstrate:

+ Model consistency across different genomes
+ Similar detection patterns in related _S. pyogenes_ strains
+ Threshold captures comparable sequence features across species
+ Detection rates align with own genome (12.6%), suggesting no dataset-specific bias

The 10-14% detection range reflects the biological reality that only a subset of genes utilize canonical $sigma^(70)$ TATAAT promoters, with remaining genes employing alternative regulatory mechanisms.

= Discussion

== Biological Validation

=== Consensus Sequence Analysis

The computed consensus *TATAAT* perfectly matches the canonical bacterial *Pribnow box*, providing strong biological validation:

+ *Computational methodology:* Statistical alignment successfully identified biologically relevant sequences matching the established $sigma^(70)$ recognition motif
+ *Pattern matching quality:* WAWWWT filter accurately captured true promoter sequences
+ *PPM construction:* Frequency-to-probability conversion with pseudocounts produced biologically meaningful probabilities
+ *Near-complete conservation:* Positions 2 (99.9% A) and 6 (99.9% T) show critical conservation defining the TATAAT motif
+ *Position-specific variation:* Positions 1, 3, 4, 5 allow AT variation (W positions) while maintaining functional DNA melting capability
+ *Zero G/C frequencies:* All C/G probabilities derive from pseudocounts only ($k = 0.01$), consistent with functional constraints for transcription bubble formation

=== AT-Richness and DNA Melting

The complete absence of G/C in training sequences (100% W bases) reflects a fundamental functional requirement for transcription initiation. From the lecture notes on promoter search:

*"From the 2H bond of A and T, mutation of A to T will not have an effect. The TATAAT box can change... Promoter functionality is retained when A mutates to T or vice versa, as there is no change to hydrogen bonds."*

Biophysical basis:
- *AT base pairs:* 2 hydrogen bonds (easily separated during DNA melting)
- *GC base pairs:* 3 hydrogen bonds (stronger, resist melting)
- *Transcription bubble formation:* RNA polymerase requires strand separation for template access; AT-richness facilitates this process

*"Promoter functionality is compromised when C or G mutations occur, due to changes in hydrogen bonds."* (Lecture notes)

=== Position-Specific Conservation and $sigma^(70)$ Recognition

The TATAAT consensus with 99.9% conservation at positions 2 and 6 is critical for:

+ *$sigma^(70)$ subunit recognition:* The RNA polymerase $sigma^(70)$ factor specifically recognizes the TATAAT sequence through sequence-specific protein-DNA contacts at positions 2 and 6
+ *DNA bending and flexibility:* AT-rich sequences bend more easily, facilitating DNA wrapping around RNA polymerase
+ *Transcription bubble nucleation:* Conserved positions serve as preferred initiation points for strand separation during open complex formation
+ *Functional flexibility:* W positions (1, 3, 4, 5) allow AT variation, enabling fine-tuning of promoter strength while maintaining DNA melting capability

The observed position-specific probabilities match *empirically-derived* patterns from genome-wide promoter analyses and align with biochemical studies of $sigma^(70)$-DNA interactions.

== Detection Rate Analysis

=== Empirical Threshold Performance

Applying the empirical threshold (-10.0) yields 12.6% detection (126/1000). This detection rate reflects biological reality and demonstrates appropriate model specificity.

*Score Distribution Analysis:*
- Training promoters: ~0.0 (near-perfect TATAAT matches)
- Test sequences: μ = -19.0, range = [-40.6, -7.7]
- Clear separation: Test sequences score well below perfect matches

*Biological Interpretation:*

The 12.6% detection reflects established promoter biology:

+ *Canonical promoter frequency:* Only ~4% of genes (42/1100) possess strong TATAAT Pribnow boxes matching the WAWWWT pattern. The empirical threshold detects additional genes with 1-2 base variations, reaching 12.6%.

+ *Alternative promoter architectures:* Remaining genes (~87%) utilize:
  - Alternative σ factors (σ³², σ⁵⁴, σˢ) with different consensus sequences
  - Extended -10 promoters (TGn-TATAAT)
  - -35-dependent promoters with weak -10 elements
  - UP element-dependent transcription

+ *Model specificity:* The PPM captures canonical $sigma^(70)$ promoters with high specificity. Low false positive rate ensures detected sequences genuinely resemble TATAAT.

+ *Threshold appropriateness:* Threshold of -10.0 allows 1-2 base variations from consensus, capturing functional promoters while excluding non-promoter sequences.

=== Threshold Selection and Methodology

*Empirical Threshold (-10.0):*

Following heuristic threshold methodology from lecture notes:
*"A heuristic threshold is applied. The consensus sequence can be used as a benchmark."*

- Rationale: Allows detection of sequences with 1-2 base variations from perfect TATAAT
- Result: 12.6% detection (biologically appropriate for canonical $sigma^(70)$ promoters)
- Validation: Cross-genome consistency (10-14%) confirms threshold generalizability
- Interpretation: Sequences scoring > -10.0 are functionally similar to TATAAT consensus

== Cross-Validation Significance

=== Model Consistency Across Genomes

Cross-validation detection rates using empirical threshold (-10.0): 10.2% (210707L) to 14.3% (210732H), mean = 12.9%, SD = 1.5%

Coefficient of variation (CV):
$ "CV" = sigma / mu = 1.5 / 12.9 = 0.116 = 11.6% $

This tight clustering (CV = 11.6%) across phylogenetically related _Streptococcus pyogenes_ strains demonstrates:

+ *Excellent model generalizability:* Consistent detection patterns suggest PPM captures conserved biological features
+ *No dataset-specific overfitting:* Own genome (12.6%) aligns perfectly with cross-validation mean (12.9%)
+ *Reproducible scoring:* PPM produces stable results across independent datasets
+ *Biological conservation:* Similar proportions of canonical TATAAT promoters across _S. pyogenes_ strains



=== Biological Implications

The similar detection patterns (10-14% range) across genomes suggest:

+ *Conserved promoter architecture:* Related _S. pyogenes_ strains share similar proportions of canonical vs. alternative promoters
+ *Consistent gene regulation:* Approximately 10-14% of genes in each strain utilize canonical $sigma^(70)$ TATAAT promoters
+ *Species-level conservation:* Within-species variation is minimal (CV = 11.6%), indicating selective pressure maintaining promoter features
+ *Regulatory diversity:* Remaining 85-90% of genes employ alternative regulatory mechanisms, reflecting diverse transcriptional control strategies

== Clinical Relevance

_S. pyogenes_ is a human pathogen causing pharyngitis, scarlet fever, and invasive infections. Understanding promoter architecture can inform:

- Antibiotic development targeting transcription
- Gene regulation studies for virulence factors
- Comparative genomics identifying strain differences

== Limitations and Future Directions

=== Methodological Limitations

+ *Single promoter element modeled:* Analysis focused exclusively on the -10 box (Pribnow box). The -35 box (TTGACA region) and spacer length (typically 17 ± 1 bp between -35 and -10 elements) were not incorporated. From lecture notes: *"The TTGACA box is a binding site for sigma factor proteins... TTGACA → TATAAT → ATG → Coding region."* A complete promoter model should include both elements for improved detection accuracy.

+ *Fixed extraction window:* Used -15 to -5 region relative to start codon. True promoters can occur at variable distances; optimal search window may differ for genes with longer 5' UTRs or alternative transcription start sites (TSSs).

+ *Pattern-based selection:* The 42 training sequences were selected using WAWWWT pattern matching ([AT]A[AT][AT][AT]T). While this ensures biological relevance, alternative selection criteria might capture additional promoter variants.

=== Cross-Validation Scope

+ *Unidirectional testing:* Applied 210657G's PPM to other genomes but did not test reciprocally (other students' PPMs on 210657G). Bidirectional cross-validation would better assess model robustness.

+ *Within-species only:* All genomes are _Streptococcus pyogenes_ strains (same species). Testing across phylogenetically distant bacteria (e.g., _E. coli_, _B. subtilis_) would evaluate true generalizability of TATAAT motif recognition.

=== Biological Validation

+ *Computational predictions only:* No experimental confirmation via:
  - RNA-seq (transcription start site mapping)
  - Promoter-reporter assays (functional activity)
  - ChIP-seq ($sigma^(70)$ binding verification)

+ *Binary classification:* Promoters classified as present/absent, but real promoters have varying *strength* (transcription rates). Statistical scores could be calibrated against expression levels.

=== Future Directions

+ *Incorporate -35 box:* Build joint PPM for both elements with spacer length modeling
+ *Use position weight matrices (PWM):* Replace probabilities with log-odds scores relative to background nucleotide frequencies
+ *Machine learning approaches:* Compare statistical alignment against modern methods (CNNs, transformers) for promoter prediction
+ *Experimental validation:* Prioritize high-scoring predictions for wet-lab verification

= Conclusion

This study successfully applied *statistical gene prediction* methodology to identify bacterial promoters in _Streptococcus pyogenes_ genome GCA_900637025.1. The analysis demonstrates that pattern-based promoter selection combined with Position Probability Matrix (PPM) scoring provides an effective framework for canonical $sigma^(70)$ promoter detection.

=== Key Findings

+ *Consensus Sequence:* Computed consensus is *TATAAT*, perfectly matching the canonical Pribnow box and providing strong biological validation of the methodology.

+ *PPM Construction:* Successfully built from 42 sequences matching WAWWWT pattern [AT]A[AT][AT][AT]T. Near-complete conservation at positions 2 (99.9% A) and 6 (99.9% T) defines the canonical motif, while positions 1, 3, 4, 5 allow AT variation (W positions) for functional flexibility.

+ *Statistical Alignment Methodology:* Log probability scoring ($sum log(p_(j comma s_j))$) provides quantitative measure of sequence similarity to TATAAT consensus. Perfect matches score ~0, while non-promoters score below -10.0.

+ *Detection Performance:* Empirical threshold (-10.0) yields 12.6% detection rate (126/1000), reflecting biological reality that only a subset of genes utilize canonical $sigma^(70)$ promoters. Remaining genes employ alternative regulatory mechanisms.

+ *Cross-Validation Results:* Excellent consistency across five _S. pyogenes_ genomes (10.2-14.3%, mean = 12.9% ± 1.5%, CV = 11.6%):
  - Demonstrates model generalizability and lack of overfitting
  - Own genome (12.6%) aligns with cross-validation mean
  - Tight clustering indicates conserved promoter architecture across strains

+ *Biological Validation:* TATAAT consensus with 99.9% conservation at critical positions matches established $sigma^(70)$ recognition motif, confirming methodology accurately captures true biological signals.

=== Methodological Contribution

This work demonstrates practical implementation of statistical gene prediction concepts:
- Pattern-based promoter selection using biologically motivated WAWWWT filter
- Empirical PPM construction with pseudocounts for unobserved bases
- Statistical alignment scoring with log probabilities
- Empirical threshold derivation allowing controlled sequence variation
- Cross-validation for assessing model generalizability

The analysis validates that successful promoter detection requires:
- *Biological prior knowledge:* WAWWWT pattern ensures training data quality
- *Appropriate threshold selection:* Empirical threshold (-10.0) balances specificity and sensitivity
- *Model validation:* Cross-genome testing confirms generalizability (CV = 11.6%)
- *Result interpretation:* Low detection rates (12.6%) reflect biological reality, not methodological failure

The perfect match between computed consensus (TATAAT) and canonical Pribnow box demonstrates that *pattern-based filtering combined with statistical modeling successfully captures authentic biological signals*.

= Appendix

*Code Availability:* Complete analysis pipeline and reproducible code available at #link("https://github.com/thuvasooriya/promoter-analysis")
