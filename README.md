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


**Programs**:

| Program | Tested version | Download / Release page |
|---------|---------------|-------------------------|
| Python | **3.11.5** | https://www.python.org/downloads/release/python-3115/ |
| R | **4.5.0** | https://cran.r-project.org/bin/ |
| IQ‑TREE 2 | **2.4.0** | https://iqtree.github.io/ |
| PRANK | **v.170427** | https://github.com/ariloytynoja/prank-msa |
| HyPhy | **2.5.62** | https://github.com/veg/hyphy |
| WebLogo 3 | **3.9.0** | https://github.com/gecrooks/weblogo |
| IUPred3 | **3** | https://iupred3.elte.hu/ |
| Jalview | **2.11.4.0** | https://www.jalview.org/download/ |

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
## Installation & Environment

<details>

1. Clone the repository:

   ```bash
   git clone https://github.com/thp42/Phylogenicity.git
   cd Phylogenicity
   ```

2. Install all dependencies, libraries and packages

3. Confirm that external binaries such as ```iqtree2```, ```prank```, and ```hyphy``` execute from the command line.

</details>


---

## Features 

<details>

   - **Orthologue retrieval** – fetch your protein family from OrthoDB.
   - **High‑quality MSA** – automated alignment with DECIPHER::AlignSeqs.
   - **Maximum‑likelihood trees** – builds robust WAG/GTR trees in IQ‑TREE 2 with 1000 bootstraps.
   - **Isoform recognition** – filters, renames, and colour‑codes isoforms for downstream plotting.
   - **Tree statistics & pruning** – suggests branch‑length cut‑offs, cleans noisy tips, and recalculates support.
   - **Motif discovery** – RegEx + IUPRED/ANCHOR/PSIPRED filters, followed by iterative PSSM scoring (PAM30 / BLOSUM62).
   - **Taxonomic annotation** – determines MRCA of motif‑bearing sequences and visualises gains/losses on a taxon‑based tree.
   - **Sequence logos** – generates WebLogo PNGs for each motif class.
   - **Codon‑level analysis** – codon alignment with Prank, dN/dS with HyPhy FEL, and ancestral sequence reconstruction.

</details>

---

## Instructions 
Follow the **“Methods Phylogenicity.docx”** stored in this repository. 

---





## Citation 

If you use these python scripts in published research, please cite:
- Manuscript
  - **Discovery of a new evolutionarily conserved short linear F-actin binding motif**: Themistoklis Paraschiakos, Biao Yuan, Kostiantyn Sopelniak, Michael Bucher, Lisa Simon, Ksenija Zonjic, Dominic Eggers, Franziska Selle, Jing Li, Stefan Linder, Thomas C. Marlovits, Sabine Windhorst. bioRxiv 2025.04.16.649135; [doi: https://doi.org/10.1101/2025.04.16.649135](https://www.biorxiv.org/content/10.1101/2025.04.16.649135v1)
- Tools:
  -   **DECIPHER**: Wright ES. DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics. 2015 Oct 6;16:322. doi: 10.1186/s12859-015-0749-z.
  -   **IQ-TREE 2**: Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol Biol Evol. 2020 May 1;37(5):1530-1534. doi: 10.1093/molbev/msaa015. Erratum in: Mol Biol Evol. 2020 Aug 1;37(8):2461. doi: 10.1093/molbev/msaa131.
  -   **PRANK**: Löytynoja A. Phylogeny-Aware Alignment with PRANK and PAGAN. Methods Mol Biol. 2021;2231:17-37. doi: 10.1007/978-1-0716-1036-7_2.
  -   **HyPhy**: Kosakovsky Pond SL, Poon AFY, Velazquez R, Weaver S, Hepler NL, Murrell B, Shank SD, Magalis BR, Bouvier D, Nekrutenko A, Wisotsky S, Spielman SJ, Frost SDW, Muse SV. HyPhy 2.5-A Customizable Platform for Evolutionary Hypothesis Testing Using Phylogenies. Mol Biol Evol. 2020 Jan 1;37(1):295-299. doi: 10.1093/molbev/msz197.
  -   **ETE3**: Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol Biol Evol. 2016 Jun;33(6):1635-8. doi: 10.1093/molbev/msw046.
  -   **WebLogo 3**: Crooks GE, Hon G, Chandonia JM, Brenner SE. WebLogo: a sequence logo generator. Genome Res. 2004 Jun;14(6):1188-90. doi: 10.1101/gr.849004. 
  -   **Jalview**: Waterhouse AM, Procter JB, Martin DM, Clamp M, Barton GJ. Jalview Version 2--a multiple sequence alignment editor and analysis workbench. Bioinformatics. 2009 May 1;25(9):1189-91. doi: 10.1093/bioinformatics/btp033.
  -   **PSSMSearch**: Krystkowiak, I., Manguy, J., & Davey, N. E. (2018). PSSMSearch: A server for modeling, visualization, proteome-wide discovery and annotation of protein motif specificity determinants. Nucleic Acids Research, 46(W1), W235–W241. https://doi.org/10.1093/nar/gky426
  -   **ANCHOR**: Dosztányi, Z., Mészáros, B., & Simon, I. (2009). ANCHOR: Web server for predicting protein binding regions in disordered proteins. Bioinformatics (Oxford, England), 25(20), 2745–2746. https://doi.org/10.1093/bioinformatics/btp518
  -   **IUPRED**:Erdős, G., Pajkos, M., & Dosztányi, Z. (2021). IUPred3: Prediction of protein disorder enhanced with unambiguous experimental annotation and visualization of evolutionary conservation. Nucleic Acids Research, 49(W1), W297–W303. https://doi.org/10.1093/nar/gkab408

---
