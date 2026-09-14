# Structural variants shape recurrent chromatin contact divergence in humans

This repository contains the custom computational code used in the study:

**Structural variants shape recurrent chromatin contact divergence in humans**

Lingbin Ni *et al.*

> **Manuscript status:** in preparation / under review
> **Repository:** `pangenome_chromatin_associated_SVs`

---

## Overview

Genetic variation can alter the linear organization of the genome and, consequently, the interpretation of chromatin contacts when sequencing data are analyzed against a single reference genome. In this study, we use donor-specific, haplotype-resolved genome assemblies together with population-scale chromatin contact maps to characterize structural variation associated with recurrent changes in three-dimensional genome organization.

## Repository organization

The code is organized according to the major analytical stages of the study. Rather than representing a single end-to-end software pipeline, each directory contains scripts used for a specific component of the analyses described in the manuscript.

```text
pangenome_chromatin_associated_SVs/
│
├── hic_processing/
├── contact_matrix/
│
├── mapping_difference_analysis/
├── mapping_shift_analysis/
├── mapping_structure_analysis/
├── chromatin_structure/
│
├── contact_difference_hotspot_identification/
├── contact_difference_hotspot_contact/
├── contact_difference_hotspot_analysis/
│
├── structure_associated_SV_panCDH/
├── structure_associated_SV_identification/
├── structure_associated_SV_feature/
│
├── expression_analysis_mapping/
├── expression_analysis_integration/
├── expression_analysis_pattern/
├── expression_analysis_correlation/
│
├── evolution_analysis_mapping/
├── evolution_analysis_feature/
│
└── README.md
```

For clarity, these directories can be grouped into five major analysis stages corresponding to the progression of the manuscript.

---

## Analysis modules

### 1. Hi-C processing and reference-dependent mapping analyses

These scripts process chromatin contact data and evaluate how the choice of genome representation influences contact mapping and inferred chromatin organization.

| Directory                      | Description                                                                                                    |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------- |
| `hic_processing/`              | Processing of Hi-C/Omni-C sequencing data and generation of filtered chromatin contact pairs.                  |
| `contact_matrix/`              | Generation and processing of chromatin contact matrices used for downstream analyses.                          |
| `mapping_difference_analysis/` | Comparison of chromatin contacts obtained from donor-specific assemblies and the shared reference genome.      |
| `mapping_shift_analysis/`      | Quantification of changes in inferred genomic positions and contact distances between mapping contexts.        |
| `mapping_structure_analysis/`  | Evaluation of mapping-dependent differences in higher-order chromatin organization.                            |
| `chromatin_structure/`         | Identification and analysis of chromatin structural features, including domain boundaries and chromatin loops. |

Together, these analyses evaluate the extent to which a shared linear reference preserves chromatin contact identity while altering the genomic coordinates or inferred separation of interacting loci.

---

### 2. Identification and characterization of contact divergence hotspots

These scripts identify genomic regions showing recurrent excesses of haplotype-specific chromatin contacts and characterize their properties across the population.

| Directory                                    | Description                                                                                                      |
| -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `contact_difference_hotspot_identification/` | Statistical identification of contact divergence hotspots (CDHs) from haplotype-resolved chromatin contact maps. |
| `contact_difference_hotspot_contact/`        | Contact-level analyses of shared and haplotype-specific interactions associated with CDHs.                       |
| `contact_difference_hotspot_analysis/`       | Population-level characterization, recurrence analysis, and genomic annotation of CDHs.                          |

CDHs are defined as genomic regions with an excess of haplotype-specific chromatin contacts relative to the local background expectation. These analyses establish the population-scale distribution and recurrence of chromatin contact divergence.

---

### 3. Identification of chromatin-structure-associated structural variants

These scripts integrate recurrent CDHs with haplotype-resolved structural variation to identify structural variants associated with chromatin contact divergence.

| Directory                                 | Description                                                                                 |
| ----------------------------------------- | ------------------------------------------------------------------------------------------- |
| `structure_associated_SV_panCDH/`         | Integration of structural variants with recurrent population-level CDH regions.             |
| `structure_associated_SV_identification/` | Statistical identification of chromatin-structure-associated structural variants (CSA-SVs). |
| `structure_associated_SV_feature/`        | Genomic and regulatory characterization of identified CSA-SVs.                              |

These analyses prioritize a subset of structural variants whose genotypes are associated with recurrent chromatin contact divergence across individuals.

---

### 4. Integration with gene-expression variation

These scripts integrate CSA-SV genotypes and chromatin organization with matched transcriptomic data.

| Directory                          | Description                                                                                     |
| ---------------------------------- | ----------------------------------------------------------------------------------------------- |
| `expression_analysis_mapping/`     | Processing and mapping of matched transcriptomic data for haplotype-aware expression analyses.  |
| `expression_analysis_integration/` | Integration of CSA-SV genotypes with gene-expression measurements and cis-association analyses. |
| `expression_analysis_pattern/`     | Characterization of genotype-dependent and allele-dosage-associated expression patterns.        |
| `expression_analysis_correlation/` | Correlation analyses between local chromatin organization and gene expression.                  |

These analyses evaluate whether structural variants associated with recurrent chromatin contact divergence are also associated with variation in nearby gene expression.

---

### 5. Evolutionary analysis of CSA-SVs

These scripts use comparative genomic information to investigate the evolutionary origin and genomic properties of CSA-SVs.

| Directory                     | Description                                                                                                |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `evolution_analysis_mapping/` | Comparative mapping and inference of structural-variant states across human and non-human primate genomes. |
| `evolution_analysis_feature/` | Characterization of the evolutionary origins and genomic features of CSA-SVs.                              |

These analyses distinguish human-derived from ancestral structural variants and evaluate the evolutionary properties of structural variants associated with chromatin contact divergence.

---

## Relationship to manuscript analyses

The code organization broadly follows the order of the major analyses presented in the manuscript:

| Manuscript analysis                                           | Relevant directories                                                                                                                                   |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Haplotype-resolved mapping and reference-dependent distortion | `hic_processing/`, `contact_matrix/`, `mapping_difference_analysis/`, `mapping_shift_analysis/`, `mapping_structure_analysis/`, `chromatin_structure/` |
| Identification of recurrent chromatin contact divergence      | `contact_difference_hotspot_identification/`, `contact_difference_hotspot_contact/`, `contact_difference_hotspot_analysis/`                            |
| Identification and characterization of CSA-SVs                | `structure_associated_SV_panCDH/`, `structure_associated_SV_identification/`, `structure_associated_SV_feature/`                                       |
| CSA-SVs and gene-expression variation                         | `expression_analysis_mapping/`, `expression_analysis_integration/`, `expression_analysis_pattern/`, `expression_analysis_correlation/`                 |
| Evolutionary analysis of CSA-SVs                              | `evolution_analysis_mapping/`, `evolution_analysis_feature/`                                                                                           |

Individual scripts within each directory correspond to specific analyses, summary statistics, and visualizations described in the Methods and figure legends of the manuscript.

---

## Scope of the repository

This repository contains the custom analysis code developed for this study. It is intended to document and reproduce the principal computational analyses rather than provide a standalone software package.

The complete cohort-scale analyses require access to the corresponding genomic datasets, donor-specific assemblies, chromatin contact data, transcriptomic data, and substantial computational resources. Large genomic input files and intermediate results are therefore not distributed through this repository.

## Versioning

The `main` branch may contain updates made during manuscript preparation and revision.

The exact code version associated with the published study will be preserved as a tagged GitHub release, for example:

```text
v1.0.0
```

A permanent archival copy with a DOI will be provided upon publication.

---

## Citation

If you use code or methods from this repository, please cite:

> Ni L. *et al.* **Structural variants shape recurrent chromatin contact divergence in humans.**
> <Journal information / DOI to be added upon publication>

Citation metadata will also be provided through:

```text
CITATION.cff
```

---

## License

This repository is distributed under the terms described in the `LICENSE` file.

Third-party software and external datasets remain subject to their respective licenses and terms of use.

---

## Contact

For questions regarding the code or analyses, please open a GitHub issue or contact:

**Lingbin Ni**
University of Washington
Eichler Laboratory

---

## Disclaimer

This repository contains research code developed for the analyses described in the accompanying manuscript. It is provided primarily to support transparency and reproducibility of the published work and should not be interpreted as a production software package.
