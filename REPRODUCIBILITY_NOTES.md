# Reproducibility notes

This repository contains analysis code rather than primary sequencing data, controlled genotype data or large processed objects. Obtain the corresponding inputs through the data repositories and access routes cited in the manuscript, then set the input paths described in each script and in [environment/README.md](environment/README.md).

Some modules consume processed matrices, metadata, model objects or result tables produced by earlier modules. SCENIC also needs cisTarget resources; QTL analysis needs separately installed OmiGA, GTEx/qtl and CattleCell-GTEx helper code. These dependencies and their input schemas should be preserved when using the scripts.

The scripts are independent analysis modules, not a single automated pipeline. Static validation checks syntax and repository-internal references; the complete analyses have not been run on the manuscript data within this code repository.
