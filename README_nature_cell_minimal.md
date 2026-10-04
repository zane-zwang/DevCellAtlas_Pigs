# DevCellAtlas_Pigs

Code accompanying a multi-tissue developmental single-nucleus transcriptomic atlas of the pig (*Sus scrofa*).

The study profiles transcriptional development across six developmental stages (E55, E90, P0, P30, P90 and P180) and seven tissues (adipose tissue, cerebrum, duodenum, heart, hypothalamus, liver and skeletal muscle), integrating cell-atlas construction, developmental transcriptional dynamics, immune development, developmental clocks, genetic regulation, transcriptome age and cross-species cortical analyses.

## Repository

```text
DevCellAtlas_Pigs/
├── scripts/
│   ├── 01_atlas_annotation/
│   ├── 02_developmental_dynamics/
│   ├── 03_immune/
│   ├── 04_developmental_clock/
│   ├── 05_bulk_eqtl/
│   ├── 06_transcriptome_age/
│   └── 07_cross_species/
├── figures/
├── environment/
├── docs/
└── REPRODUCIBILITY_NOTES.md
```

| Directory | Analysis |
|---|---|
| `scripts/01_atlas_annotation/` | Atlas construction, integration and cell-type annotation |
| `scripts/02_developmental_dynamics/` | Developmental differential expression and regulon analysis |
| `scripts/03_immune/` | Immune-cell states, tissue specialization and differential abundance |
| `scripts/04_developmental_clock/` | Cell-type-specific and bulk developmental clocks |
| `scripts/05_bulk_eqtl/` | Developmental bulk and cell-type-resolved eQTL analyses |
| `scripts/06_transcriptome_age/` | Transcriptome age index analyses |
| `scripts/07_cross_species/` | Cross-species cortical alignment and OPC–OL trajectory analyses |
| `figures/` | Main and supplementary figure-generation scripts |

Scripts are numbered according to their logical or execution order where dependencies exist.

## Usage

This repository contains analysis modules rather than a single end-to-end workflow.

Set the repository and analysis-data locations before running the relevant module:

```bash
export DEVPIGATLAS_ROOT=/path/to/DevCellAtlas_Pigs
export DEVPIGATLAS_ANALYSIS_DIR=/path/to/analysis_data
```

Module-specific dependencies and path requirements are described in [`environment/README.md`](environment/README.md).

## Data availability

Raw sequencing data, processed datasets and large intermediate analysis objects are not hosted in this repository. Accession information is provided in the accompanying manuscript and associated public data repositories.

## Software

Analyses were performed using R, Python and Bash, with major components including Seurat, Scanpy, scVI-tools, pySCENIC, Milo, edgeR, glmnet, OmiGA, mashR, Bisque/bMIND, myTAI, CAMEX, Slingshot, GenEra and DIAMOND.

See [`environment/DEPENDENCIES.tsv`](environment/DEPENDENCIES.tsv) for the software inventory.

## Reproducibility

The repository provides the analysis and figure-generation code used for the manuscript. Large intermediate objects and external reference resources are distributed separately or must be generated from the corresponding upstream analyses.

Additional notes are provided in [`REPRODUCIBILITY_NOTES.md`](REPRODUCIBILITY_NOTES.md).

## Citation

If you use this code, please cite the accompanying manuscript.

Citation details will be updated upon publication.
