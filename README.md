# Phylogenicity – Phylogenetic & Motif‑Evolution Toolkit

Phylogenicity is a collection of R, Python, and shell helpers that automate the
evolutionary analysis of short‑linear motifs (SLiMs) across large protein
families.  Starting from raw orthologue FASTA files, the pipeline produces
multiple‑sequence alignments, maximum‑likelihood trees, isoform‑specific
clustering, motif extraction and scoring, MRCA determination, sequence logos,
and codon‑based dN/dS tests.

Detailed step‑by‑step instructions (including all 25 code blocks) are provided
in **“Methods Phylogenicity.docx”** stored in this repository.  The README gives
only a high‑level overview; for exact commands, parameters, and expected
outputs **please consult that Word file**.


---

## System Requirements

<details>

---

### Software

#### **OS Requirements**
This package is supported for Linux. The package has been tested on the following systems:

  - Linux: Ubuntu 20.04

#### **Dependencies**
The following packages/programs should be installed. The typical installation time should take some minutes. 

| Component | Minimum Version | Purpose |
|-----------|-----------------|---------|
| **Python** | 3.11 | motif extraction, PSSM, tree parsing |
| **R** | 4.2 | sequence alignment with **DECIPHER** |
| **IQ‑TREE 2** | 2.2 | maximum‑likelihood tree building |
| **Prank** | 170427 | codon‑aware alignment |
| **HyPhy** | 2.5 | dN/dS & ancestral reconstruction |
| **WebLogo 3** | 3.8 | sequence‑logo rendering |

**Python packages**:
```
numpy v.1.26.2
pandas v.2.2.2
biopython v.1.83
seaborn v.0.13.1
matplotlib v.3.9.2
ete3 v.3.1.2
```

**R packages** (install from CRAN/Bioconductor):
```
DECIPHER (v.2.28)
```
---

### Hardware Requirements

#### **Recommended System**

| Component | Specification | Notes |
|-----------|---------------|-------|
| **CPU** | ≥ 6 cores / 12 threads | Needed for Sequence Alignment. |
| **RAM** | ≥ 32 GB (128 GB preferred) | Needed for Sequence Alignment. |
| **GPU** | None | No GPU needed here. |
| **Storage** | ≥ 1 TB SSD (NVMe recommended) | Databases and intermediates benefit from fast I/O. |


#### **Database Storage Requirements**

| Database                       | Approximate Size |
|--------------------------------|------------------|
| Files (MSA, logs, Analysis)    | Variable (~1 GB per protein) |


</details>

---
