<div align="center">

# DevCellAtlas_Pigs

### A multi-tissue developmental single-nucleus transcriptomic atlas of the pig

*Analysis code accompanying a developmental atlas of* ***Sus scrofa*** *across tissues, cell types and developmental stages.*

[Overview](#overview) · [Analysis](#analysis-overview) · [Software](#software) · [Data](#data-availability) · [Citation](#citation-and-license)

</div>

---

## Overview

This repository contains the analysis and figure-generation code for a multi-tissue developmental single-nucleus RNA-sequencing study of the pig (*Sus scrofa*), spanning embryonic to postnatal development. The study investigates developmental transcriptional dynamics, immune maturation, transcriptomic clocks, genetic regulation, transcriptome age and cross-species cortical development.

<div align="center">

**7 tissues · 6 developmental stages · Single-nucleus RNA sequencing**

Adipose · Cerebrum · Duodenum · Heart · Hypothalamus · Liver · Skeletal muscle

**Embryonic** &nbsp; E55 · E90 &nbsp;&nbsp; | &nbsp;&nbsp; **Postnatal** &nbsp; P0 · P30 · P90 · P180

</div>

The repository serves as a companion to the manuscript and an archive of the associated analysis and figure-generation code, rather than a standalone executable pipeline.

## Analysis overview

The developmental atlas provides a common reference for the complementary research modules below. Each linked title leads to the corresponding analysis code; figure numbers refer to the manuscript. Additional bulk, genetic or cross-species inputs are used where appropriate.

<div align="center">
<table align="center" width="100%">
<tr>
<td align="center" valign="middle" colspan="2">
<sub>FIGURE 1</sub><br>
<strong><a href="scripts/01_atlas_annotation/">Atlas construction and annotation</a></strong><br>
<sub>QC · AnnData · scVI · Graph construction · Annotation</sub>
</td>
</tr>
<tr><td colspan="2" align="center"><sub><strong>COMPLEMENTARY RESEARCH MODULES</strong></sub></td></tr>
<tr>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 2</sub><br>
<strong><a href="scripts/02_developmental_dynamics/">Developmental dynamics</a></strong><br>
<sub>Differential expression · Enrichment · SCENIC · CSI</sub>
</td>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 3</sub><br>
<strong><a href="scripts/03_immune/">Immune development</a></strong><br>
<sub>Immune integration · Subtypes · Maturation · Milo</sub>
</td>
</tr>
<tr>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 4</sub><br>
<strong><a href="scripts/04_developmental_clock/">Developmental clocks</a></strong><br>
<sub>Cell-type-specific clocks · Bulk validation</sub>
</td>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 5</sub><br>
<strong><a href="scripts/05_bulk_eqtl/">Developmental genetic regulation</a></strong><br>
<sub>Bulk / interaction / cell-type eQTL · mashR · Colocalization</sub>
</td>
</tr>
<tr>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 6</sub><br>
<strong><a href="scripts/06_transcriptome_age/">Transcriptome age</a></strong><br>
<sub>Gene age · Pseudobulk · TAI</sub>
</td>
<td align="center" valign="middle" width="50%">
<sub>FIGURE 7</sub><br>
<strong><a href="scripts/07_cross_species/">Cross-species cortical development</a></strong><br>
<sub>Integration · DTW · OPC–OL · Human-centered comparisons</sub>
</td>
</tr>
</table>
</div>

Scripts are numbered according to their intended logical or execution order where dependencies exist.

**Figure generation.** Scripts for the main and supplementary figures are provided under [`figures/`](figures/). Final multi-panel assembly and graphical layout were performed separately.

## Software

Analyses were implemented in **R**, **Python** and **Bash**. Major software components include:

**Single-cell genomics and regulatory analysis** — Seurat, Scanpy, scVI-tools, edgeR, pySCENIC, Milo, Slingshot.  
**Statistical modeling and genetics** — glmnet, OmiGA, mashR, Bisque/bMIND.  
**Evolutionary and comparative genomics** — GenEra, DIAMOND, myTAI, CAMEX.

Software requirements vary by analysis module. Refer to the corresponding scripts and the manuscript Methods for implementation details and reported versions.

## Data availability

Raw sequencing data, processed datasets and large intermediate analysis objects are not hosted in this repository. Data accession information is provided in the accompanying manuscript and associated public repositories.

Some analyses additionally require external reference resources, genome annotations, cisTarget databases or third-party QTL utilities.

## Citation and license

If you use this repository or adapt the accompanying analysis code, please cite the associated manuscript.

> Citation details will be added upon publication.

See [LICENSE](LICENSE) for the repository's licensing terms.
