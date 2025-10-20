# BM4322 Genomic Signal Processing

**2021 Batch - Academic Year 2025/2026 - Semester 7**

**Assignment (Individual) - 50% of Final Grade**

---

According to Liu et al. 2011[^1], the σ² subunit of bacterial RNA polymerase is capable of detecting promoters of the form **WAWWWT** where **W** represents a **T** or **A**. These promoters typically occur **10 bases upstream**. Based on this postulate perform the following computational operations on the sequenced DNA of the given organism.

## Tasks

### 1. Position Probability Matrix Construction

Locate **1100 genes** of the organism that have been predicted based upon homology. Select the region from **15 to 5 bases upstream** for each gene. Using **100** of these selected regions manually extract the **6 bases** that most likely form the promoter. If the region does not have at least **6 consecutive Ws** reject it. Construct the **position probability matrix (PPM)** of the promoter from this data taking suitable heuristic probability values for bases **C** and **G**.

### 2. Statistical Alignment on Remaining Regions

For the remaining **1000 regions** perform a statistical alignment using the PPM of (1) and determine the presence or absence of promoters within them.

### 3. Cross-Validation with Other Students' Data

Obtain the **1000 samples** of every other student of the class and repeat (2) with the PPM of (1).

---

## Genome Assignments

The genome assigned to each student is as follows:

| Student | Genome          |
| ------- | --------------- |
| 210079K | GCA_001457635.1 |
| 210179R | GCA_019048645.1 |
| 210504L | GCA_900636475.1 |
| 210657G | GCA_900637025.1 |
| 210707L | GCA_900475505.1 |
| 210732H | GCA_019046945.1 |

---

## Important Information

**Due Date:** 2025.10.17

**Weight:** This assignment accounts for **50% of the module assessment**

**Instructor:** Upeka Premaratne  
**Email:** upeka@uom.lk  
**Contact:** 0719538433 (voice and WhatsApp)

**Assignment Released:** 2025.09.19

---

[^1]: Liu, X., Bushnell, D. A., & Kornberg, R. D. (2011). Lock and key to transcription: σ-DNA interaction. _Cell_, _147_(6), 1218-1219

