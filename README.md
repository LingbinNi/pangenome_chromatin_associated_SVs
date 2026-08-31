# Structural variants shape recurrent chromatin contact divergence in humans

This repository contains the custom computational code used in the study:

**Structural variants shape recurrent chromatin contact divergence in humans**

Lingbin Ni *et al.*

> **Manuscript status:** in preparation / under review
> **Repository:** `pangenome_chromatin_associated_SVs`

---

## Overview

Genetic variation can alter the linear organization of the genome and, consequently, the interpretation of chromatin contacts when sequencing data are analyzed against a single reference genome. In this study, we use donor-specific, haplotype-resolved genome assemblies together with population-scale chromatin contact maps to characterize structural variation associated with recurrent changes in three-dimensional genome organization.

The analyses implemented in this repository include:

* haplotype-resolved mapping of chromatin contact data to donor-specific assemblies;
* comparison of chromatin contacts obtained from donor-specific and reference-based mappings;
* quantification of haplotype-specific chromatin contact divergence;
* identification of **contact divergence hotspots (CDHs)**;
* association of structural variants with recurrent CDHs;
* identification of **chromatin-structure-associated structural variants (CSA-SVs)**;
* integration of CSA-SVs with chromatin features and gene-expression variation;
* evolutionary analyses of human-derived CSA-SVs;
* selected genotype–phenotype association analyses; and
* generation of figures and summary statistics reported in the manuscript.

The primary human dataset comprises haplotype-resolved genome assemblies and chromatin contact maps from the Human Pangenome Reference Consortium (HPRC).

---

## Analysis workflow

The overall analysis framework is:

```text
Chromatin contact sequencing
        │
        ▼
Donor-specific / haplotype-resolved mapping
        │
        ▼
Contact filtering and coordinate comparison
        │
        ▼
Haplotype-specific contact quantification
        │
        ▼
Contact divergence hotspot (CDH) detection
        │
        ▼
Structural variant association
        │
        ▼
Chromatin-structure-associated SVs (CSA-SVs)
        │
        ├── Chromatin architecture
        ├── Gene expression
        ├── Evolutionary analyses
        └── Phenotypic associations
```

---

## Repository structure

```text
pangenome_chromatin_associated_SVs/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
│
├── config/
│   └── example_config.yaml
│
├── scripts/
│   │
│   ├── 01_mapping/
│   │   └── Scripts for reference- and donor-specific mapping
│   │
│   ├── 02_contact_comparison/
│   │   └── Comparison of contact coordinates and mapping contexts
│   │
│   ├── 03_haplotype_contacts/
│   │   └── Identification and quantification of haplotype-specific contacts
│   │
│   ├── 04_CDH_detection/
│   │   └── Statistical detection and characterization of CDHs
│   │
│   ├── 05_SV_association/
│   │   └── Association analyses between structural variants and CDHs
│   │
│   ├── 06_expression/
│   │   └── Integration of CSA-SV genotypes with gene expression
│   │
│   ├── 07_evolution/
│   │   └── Evolutionary and comparative genomic analyses
│   │
│   ├── 08_PheWAS/
│   │   └── Selected phenotype association analyses
│   │
│   └── 09_figures/
│       └── Scripts used to generate manuscript figures
│
├── example/
│   ├── input/
│   └── expected_output/
│
└── docs/
    └── workflow.md
```

The exact directory organization may evolve during manuscript revision. The version associated with the published manuscript will be archived as a fixed release.

---

## Major analysis modules

### 1. Donor-specific chromatin contact mapping

Chromatin contact sequencing reads are mapped independently to the two haplotype-resolved assemblies of each donor and, where appropriate, to the T2T-CHM13 reference genome.

This module contains scripts for:

* read mapping;
* contact-pair parsing and filtering;
* haplotype-specific contact assignment;
* generation of contact matrices; and
* comparison between donor-specific and reference-based mappings.

Relevant directory:

```text
scripts/01_mapping/
```

---

### 2. Reference-dependent contact differences

Contacts detected using donor-specific assemblies are compared with those obtained after mapping the same sequencing data to a common reference genome.

These analyses quantify:

* recovery of valid contact pairs;
* changes in inferred genomic distance;
* local mapping-dependent contact differences; and
* differences in loop calls between mapping contexts.

Relevant directory:

```text
scripts/02_contact_comparison/
```

---

### 3. Haplotype-specific chromatin contacts

Haplotype-resolved contact maps are used to quantify chromatin interactions preferentially supported by one haplotype.

Relevant directory:

```text
scripts/03_haplotype_contacts/
```

This module includes procedures for:

* identifying shared and haplotype-specific contacts;
* filtering contacts by mapping quality and genomic distance;
* summarizing haplotype-specific contact density; and
* evaluating potential mapping ambiguity.

---

### 4. Contact divergence hotspots

**Contact divergence hotspots (CDHs)** are genomic regions showing an excess of haplotype-specific chromatin contacts relative to the local background expectation.

CDH detection accounts for the relationship between shared and haplotype-specific contact coverage and for coverage-dependent variability.

Relevant directory:

```text
scripts/04_CDH_detection/
```

The statistical framework includes:

* coverage-aware modeling of expected haplotype-specific contact counts;
* estimation of coverage-dependent dispersion;
* upper-tail statistical testing;
* multiple-testing correction; and
* merging and characterization of significant windows.

See the manuscript Methods for the complete statistical definition and filtering criteria.

---

### 5. Structural variant association and CSA-SVs

Structural variants overlapping or occurring near recurrent CDHs are tested for association with chromatin contact divergence across individuals and haplotypes.

Variants satisfying the predefined statistical and recurrence criteria are referred to as:

> **chromatin-structure-associated structural variants (CSA-SVs)**

Relevant directory:

```text
scripts/05_SV_association/
```

This module includes:

* SV–CDH overlap analyses;
* enrichment analyses;
* genotype–chromatin association testing;
* permutation procedures;
* recurrence filtering; and
* CSA-SV annotation.

---

### 6. Gene-expression integration

CSA-SV genotypes are integrated with matched transcriptomic data to identify loci associated with gene-expression differences.

Relevant directory:

```text
scripts/06_expression/
```

Analyses include:

* genotype–expression association testing;
* expression dosage-pattern classification;
* chromatin–expression correlation analyses; and
* locus-level visualization.

---

### 7. Evolutionary analyses

Comparative genomic analyses are used to infer the evolutionary origin of selected CSA-SVs and evaluate chromatin divergence at orthologous loci across humans and non-human primates.

Relevant directory:

```text
scripts/07_evolution/
```

---

### 8. Phenotype association analyses

Selected CSA-SVs were evaluated for associations with human phenotypes using population-scale genotype and phenotype data.

Relevant directory:

```text
scripts/08_PheWAS/
```

Because some analyses use controlled-access datasets, individual-level genotype and phenotype data are **not distributed through this repository**.

---

## Software and computational environment

Analyses were performed using a combination of Python, R, shell scripts, and established genomics software.

Major dependencies include:

```text
Python
R
BWA
SAMtools
BCFtools
BEDTools
pairtools
cooler
cooltools
HiCExplorer
PLINK2
```

Exact software versions used for the manuscript analyses are documented in:

```text
environment.yml
```

and/or within the corresponding analysis directories.

To create the Conda environment:

```bash
conda env create -f environment.yml
conda activate pangenome-chromatin-sv
```

Some analyses require additional software or high-performance computing resources that are not distributed through Conda. These dependencies are described in the relevant scripts or documentation.

---

## Usage

Most scripts are designed to operate on preprocessed genomic files and accept input/output paths as command-line arguments.

A typical analysis can be executed as:

```bash
python scripts/<analysis_module>/<script>.py \
    --input <input_file> \
    --output <output_file>
```

or:

```bash
Rscript scripts/<analysis_module>/<script>.R \
    <input_file> \
    <output_file>
```

Scripts requiring cohort-scale processing may additionally require configuration files describing sample names, genome assemblies, reference files, or computational resources.

Example configurations are provided in:

```text
config/
```

---

## Example data

Because the complete study involves large-scale chromatin contact datasets and hundreds of haplotype-resolved genome assemblies, the full analysis cannot be reproduced using data stored directly in this GitHub repository.

Where feasible, small example inputs and expected outputs are provided in:

```text
example/
```

These files are intended to illustrate input formats and allow key custom scripts to be tested without downloading the full study dataset.

---

## Reproducing manuscript analyses

The relationship between major manuscript analyses and code modules is summarized below.

| Manuscript analysis                     | Code                             |
| --------------------------------------- | -------------------------------- |
| Donor-specific mapping                  | `scripts/01_mapping/`            |
| Reference-dependent contact differences | `scripts/02_contact_comparison/` |
| Haplotype-specific contacts             | `scripts/03_haplotype_contacts/` |
| CDH identification                      | `scripts/04_CDH_detection/`      |
| SV enrichment and CSA-SV discovery      | `scripts/05_SV_association/`     |
| CSA-SV–expression analyses              | `scripts/06_expression/`         |
| Evolutionary analyses                   | `scripts/07_evolution/`          |
| Phenome-wide association analyses       | `scripts/08_PheWAS/`             |
| Figure generation                       | `scripts/09_figures/`            |

Additional script-to-figure mappings will be provided with the final manuscript release.

---

## Data availability

The analyses in this study integrate publicly available and controlled-access genomic datasets.

Large genomic datasets, including raw chromatin contact sequencing data, genome assemblies, contact matrices, and controlled-access genotype/phenotype data, are not stored directly in this repository.

Data-access information and accession identifiers are described in the **Data availability** section of the manuscript.

Relevant resources include:

* Human Pangenome Reference Consortium datasets;
* donor-specific haplotype-resolved genome assemblies;
* matched chromatin contact datasets;
* matched transcriptomic datasets; and
* controlled-access population-scale genotype and phenotype datasets used for selected association analyses.

Users are responsible for obtaining any required data-access approvals before reproducing analyses involving controlled-access datasets.

---

## Reproducibility

The repository is intended to provide the custom computational procedures required to reproduce the principal analyses reported in the manuscript.

For reproducibility:

1. custom analysis scripts are organized according to the scientific workflow;
2. command-line parameters and input/output formats are documented where applicable;
3. software versions are recorded;
4. small example datasets are provided when redistribution is permitted; and
5. the code corresponding to the published manuscript will be archived as a versioned release.

The complete cohort-scale analyses require substantial computational resources and access to the source genomic datasets and therefore are not intended to run directly from this repository as a single end-to-end workflow.

---

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
