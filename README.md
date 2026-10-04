# DevCellAtlas_Pigs

Analysis code accompanying a multi-tissue developmental single-nucleus RNA-seq atlas of the pig (*Sus scrofa*).

The study profiles transcriptional development across six stages (E55, E90, P0, P30, P90 and P180) and seven tissues (adipose tissue, cerebrum, duodenum, heart, hypothalamus, liver and skeletal muscle). The repository contains the analysis and visualization code used for the main and supplementary results, including cell-atlas construction, developmental transcriptional dynamics, immune-cell analysis, developmental clocks, developmental eQTL analyses, transcriptome age index analysis and cross-species cortical comparisons.

## Repository structure

| Path | Description |
|---|---|
| `scripts/01_atlas_annotation/` | Quality control, AnnData construction, scVI integration, graph construction and atlas-level clustering |
| `scripts/02_developmental_dynamics/` | Adjacent-stage differential expression, functional enrichment, SCENIC and CSI analyses |
| `scripts/03_immune/` | Immune-cell integration, subtype programs, cross-tissue maturation modeling and Milo differential abundance |
| `scripts/04_developmental_clock/` | Cell-type-specific and bulk developmental clock analyses |
| `scripts/05_bulk_eqtl/` | Bulk eQTL, mashR, interaction eQTL, cell-type-resolved eQTL, enrichment and colocalization analyses |
| `scripts/06_transcriptome_age/` | Gene-age inference, pseudobulk expression processing and transcriptome age index analyses |
| `scripts/07_cross_species/` | Cross-species integration, developmental-stage alignment, OPC–OL trajectory analysis and human-centered comparisons |
| `figures/` | Code used to generate main and supplementary figure panels |
| `environment/` | Software requirements and input/path configuration notes |
| `docs/` | Script provenance and release metadata |
| `REPRODUCIBILITY_NOTES.md` | Scope and reproducibility notes for the released code |

Scripts within each analysis directory are numbered in their intended logical or execution order where applicable.

## Analysis modules and manuscript figures

| Module | Main analyses | Manuscript figure |
|---|---|---|
| Atlas construction and annotation | snRNA-seq QC, scVI integration, clustering, annotation and variance decomposition | Figure 1 |
| Developmental transcriptional dynamics | Adjacent-stage differential expression, enrichment and regulon-module analysis | Figure 2 |
| Immune development | Immune subtyping, tissue programs, maturation models and differential abundance | Figure 3 |
| Developmental clocks | Cell-type-specific LASSO clocks and bulk skeletal-muscle validation | Figure 4 |
| Developmental genetic regulation | Bulk, interaction and cell-type-resolved eQTL analyses, mashR and colocalization | Figure 5 |
| Transcriptome age | Gene-age assignment, TAI modeling and cross-species developmental comparison | Figure 6 |
| Cross-species cortical development | CAMEX integration, DTW alignment, OPC–OL trajectories and human-centered comparisons | Figure 7 |

## Usage

This repository is organized as a collection of analysis modules rather than a single end-to-end workflow. Most scripts begin from processed count matrices, metadata, genotype data or intermediate analysis objects generated as described in the manuscript.

For repository-internal references, set:

```bash
export DEVPIGATLAS_ROOT=/path/to/DevCellAtlas_Pigs
