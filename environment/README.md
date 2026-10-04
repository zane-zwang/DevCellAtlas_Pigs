# Software and inputs

The modules use R, Python and Bash. [DEPENDENCIES.tsv](DEPENDENCIES.tsv) inventories external R/Python packages and command-line programs observed in the scripts; it is not a lockfile. Version anchors supported by manuscript/environment records are Milo 1.2.0, CAMEX 0.0.1, GenEra 1.4.2, DIAMOND 2.1.10 and myTAI 2.3.5. Other package versions are left unspecified.

Set `DEVPIGATLAS_ROOT` to the repository directory for included `source()` calls. Set `DEVPIGATLAS_ANALYSIS_DIR` where a script expects external matrices, metadata, genotype files, annotations or intermediate results. SCENIC additionally needs cisTarget feather/ranking resources and pySCENIC. Set `SCENIC_SEURAT_RDS`, `SCENIC_REFERENCE_DIR`, `SCENIC_METADATA_CSV` and `SCENIC_OUTPUT_DIR`; the shell script passes matching counts, loom and ranking paths to the included helpers. All three SCENIC helpers are under `scripts/02_developmental_dynamics/helpers/`.

The adipose cell-type clock accepts `CLOCK_INPUT_RDS`. The bMIND module accepts `BMIND_REFERENCE_RDS`, `BMIND_BISQUE_FRACTIONS_CSV`, `BMIND_GENE_IDS_CSV` and `BMIND_BULK_COUNTS_CSV`; set these to the matching manuscript input files before running that module. File-level paths do not establish sample or cell-type alignment, which requires the associated metadata.

The clock-gene overlap script accepts `CLOCK_GENE_INPUT_DIR`, `CLOCK_ATLAS_METADATA_CSV`, `CLOCK_PUBLIC_DATA_DIR`, `CLOCK_IMPORTANCE_CSV`, `HUMAN_CLOCK_DATA_DIR` and `CLOCK_OVERLAP_OUTPUT_DIR`. The cross-species TAI bulk comparison accepts `TAI_BULK_INPUT_DIR`; supplementary annotation plotting accepts `ANNOTATION_INPUT_DIR`. OPC–OL Slingshot and trajectory classification share `OPC_OL_OUTPUT_DIR` and `OPC_OL_SLINGSHOT_SEURAT_RDS`. These variables select input/output locations without changing the analysis methods.

The bulk, interaction and cell-specific QTL association filters accept `QTL_BULK_ASSOCIATION_DIR`, `QTL_INTERACTION_ASSOCIATION_DIR`, `QTL_CSE_ASSOCIATION_DIR` and their corresponding `QTL_*_FILTER_OUTPUT_DIR` variables. The bulk filter output location should also be supplied to strong-pair extraction when overriding its default.

The QTL shell modules require OmiGA, PLINK, and an external GTEx/qtl expression-preparation script authored by Francois Aguet (`GTEX_PIPELINE_DIR/qtl/src/eqtl_prepare_expression_tss_gene.py`). The mashR strong-pair extraction command calls `combine_signif_pairs_tjy.py` and `extract_pairs_tjy.py` from the third-party CattleCell-GTEx pipeline (`CATTLECELL_PIPELINE_DIR/6.Tissue and Celltype sharing/`). These recovered external files are not rebranded as project-authored code or bundled here; obtain and verify their licensing and runtime versions before use. TAI gene-age preparation invokes `ncbitax2lin`, DIAMOND and GenEra.

Supply the exact input tables and object schemas expected by each script. Numerical reproduction and runtime compatibility have not been established by static checks alone.
