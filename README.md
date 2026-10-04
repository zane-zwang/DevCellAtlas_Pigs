# DevCellAtlas_Pigs analysis code

Analysis code accompanying the DevCellAtlas_Pigs manuscript, a multi-tissue developmental single-nucleus RNA-seq atlas of pigs. The repository contains independent analysis modules and figure-panel scripts. It is not a one-command end-to-end pipeline.

| Directory | Analysis | Main figure |
|---|---|---|
| `scripts/01_atlas_annotation/` | Example-sample QC, atlas scVI integration, clustering and annotation inputs | Figure 1 |
| `scripts/02_developmental_dynamics/` | Adjacent-stage differential-expression statistics, enrichment, SCENIC and CSI | Figure 2 |
| `scripts/03_immune/` | Immune integration, transcriptional programs and differential abundance | Figure 3 |
| `scripts/04_developmental_clock/` | Cell-type and bulk developmental clocks | Figure 4 |
| `scripts/05_bulk_eqtl/` | Bulk eQTL, mashR, interaction/cell-type QTL and colocalization | Figure 5 |
| `scripts/06_transcriptome_age/` | Gene ages, CPM expression summaries and TAI | Figure 6 |
| `scripts/07_cross_species/` | CAMEX, stage alignment, nine-bin OPC–OL trajectories and human-centered comparisons | Figure 7 |
| `figures/` | Main and supplementary figure-panel code | Figures 1–7 and Supplementary Figures |

The immune analysis includes 65,026 cells. Cell-type clock eligibility uses more than 100 nuclei per cell-type-stage group and representation in at least four stages. OPC–OL trajectories use nine pseudotime bins (Early 1–3, Mid 4–6, Late 7–9). snRNA-seq TAI uses CPM and the log transformation implemented in the scripts. The adjacent-stage FindMarkers script generates differential-expression statistics; manuscript Methods specify downstream selection criteria.

Start with the script for the analysis of interest. Set `DEVPIGATLAS_ROOT` to this repository for included `source()` calls. Most modules accept `DEVPIGATLAS_ANALYSIS_DIR` for external data and results; additional file-level variables are described in [environment/README.md](environment/README.md). For SCENIC, set the Seurat input and cisTarget resources, then run `03.run_pyscenic.sh`; its outputs and the CSI analysis share `SCENIC_OUTPUT_DIR`. For OPC–OL analysis, run `11.slingshot.R` before `12.trajectory_wide_classification.R`; both use `OPC_OL_OUTPUT_DIR`. QTL scripts require matched genotype, expression and phenotype inputs, OmiGA, and separately installed attributed third-party helpers. [REPRODUCIBILITY_NOTES.md](REPRODUCIBILITY_NOTES.md) describes input boundaries.

Large primary sequencing datasets and major intermediate objects are distributed through the repositories cited in the manuscript rather than GitHub. The code begins from count/expression matrices after alignment and count generation. Final figure assembly was performed in Affinity Designer. [docs/SOURCE_PROVENANCE.tsv](docs/SOURCE_PROVENANCE.tsv) lists published scripts and checksums. Static syntax and internal-reference checks do not establish numerical reproduction.
