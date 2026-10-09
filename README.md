# Structural variants shape recurrent chromatin contact divergence in humans

This repository contains custom computational code used in the study:

**Structural variants shape recurrent chromatin contact divergence in humans**

Lingbin Ni *et al.*

> **Manuscript status:** in preparation / under review
> **Repository:** `pangenome_chromatin_associated_SVs`

---

## Overview

Genetic variation can alter the linear organization of the genome and, consequently, the interpretation of chromatin contacts when sequencing data are analyzed against a single reference genome. This study uses donor-specific, haplotype-resolved genome assemblies and population-scale chromatin contact maps to characterize structural variants associated with recurrent differences in three-dimensional genome organization.

The scripts cover human and nonhuman-primate contact mapping, contact-matrix construction, chromatin-feature calling, reference-dependent contact comparisons, contact difference hotspot (CDH) identification, chromatin-structure-associated structural variant (CSA-SV) identification, and selected expression and evolutionary analyses.

The primary human dataset comprises 177 individuals with 354 haplotype-resolved assemblies and matched chromatin contact data from the Human Pangenome Reference Consortium (HPRC). Matched long-read RNA-seq data are available for 163 individuals. Individual analyses use the samples meeting their input and filtering requirements; these cohort totals are not the sample size of every test.

This is a collection of analysis scripts, rather than a single automated workflow. Intermediate-data preparation and final multi-panel figure assembly may be required. The current repository does not include the PheWAS pipeline, ancestral-state inference, or the ancestry-adjusted and donor-clustered robustness analyses described in the manuscript.

---

## Repository structure

Scripts are stored directly in the following analysis directories.

| Directory | Scripts | Purpose |
| --- | --- | --- |
| `hic_processing/` | [map_hic_to_t2t_chm13.sh](hic_processing/map_hic_to_t2t_chm13.sh), [map_hic_to_dsa_hap1.sh](hic_processing/map_hic_to_dsa_hap1.sh), [map_hic_to_dsa_hap2.sh](hic_processing/map_hic_to_dsa_hap2.sh) | Map human chromatin contact reads to CHM13 and donor-specific haplotypes. |
| `contact_matrix/` | [build_contact_matrices_t2t_chm13.sh](contact_matrix/build_contact_matrices_t2t_chm13.sh), [build_contact_matrices_dsa_hap1.sh](contact_matrix/build_contact_matrices_dsa_hap1.sh), [build_contact_matrices_dsa_hap2.sh](contact_matrix/build_contact_matrices_dsa_hap2.sh) | Build contact matrices from mapped pairs. |
| `chromatin_structure/` | [call_tads_hicexplorer.sh](chromatin_structure/call_tads_hicexplorer.sh), [call_loops_hicexplorer.sh](chromatin_structure/call_loops_hicexplorer.sh) | Call TADs, boundaries and chromatin loops. |
| `mapping_difference_analysis/` | [compare_mapping_chm13_dsa_hap1.sh](mapping_difference_analysis/compare_mapping_chm13_dsa_hap1.sh), [compare_mapping_chm13_dsa_hap2.sh](mapping_difference_analysis/compare_mapping_chm13_dsa_hap2.sh) | Compare retained read-pair IDs and contact distances between mapping targets. |
| `mapping_shift_analysis/` | [summarize_mapping_distance_differences_hap1.sh](mapping_shift_analysis/summarize_mapping_distance_differences_hap1.sh), [summarize_mapping_distance_differences_hap2.sh](mapping_shift_analysis/summarize_mapping_distance_differences_hap2.sh) | Summarize counts of contact-distance differences above specified thresholds. |
| `mapping_structure_analysis/` | [compare_loops_liftover.py](mapping_structure_analysis/compare_loops_liftover.py) | Compare loop anchors after PAF-based coordinate liftover. |
| `contact_difference_hotspot_contact/` | [extract_haplotype_specific_pairs.sh](contact_difference_hotspot_contact/extract_haplotype_specific_pairs.sh), [remap_haplotype_specific_reads.sh](contact_difference_hotspot_contact/remap_haplotype_specific_reads.sh), [build_haplotype_mappability_tracks.sh](contact_difference_hotspot_contact/build_haplotype_mappability_tracks.sh) | Extract contacts from prepared read-ID lists, remap selected reads and build mappability tracks. |
| `contact_difference_hotspot_identification/` | [count_long_range_cis_contacts.pl](contact_difference_hotspot_identification/count_long_range_cis_contacts.pl), [call_contact_difference_hotspots.R](contact_difference_hotspot_identification/call_contact_difference_hotspots.R) | Count long-range contacts and identify CDH windows and merged regions. |
| `contact_difference_hotspot_analysis/` | [test_cdh_structure_enrichment.R](contact_difference_hotspot_analysis/test_cdh_structure_enrichment.R), [test_cdh_sv_enrichment.R](contact_difference_hotspot_analysis/test_cdh_sv_enrichment.R), [calculate_sv_in_cdh_fraction.R](contact_difference_hotspot_analysis/calculate_sv_in_cdh_fraction.R) | Evaluate CDH overlap with chromatin features and SVs, and calculate the fraction of SVs overlapping CDHs. |
| `structure_associated_SV_panCDH/` | [analyze_cdh_recurrence.R](structure_associated_SV_panCDH/analyze_cdh_recurrence.R) | Identify recurrent CDH cores from projected intervals and evaluate recurrence by circular-shift permutations. |
| `structure_associated_SV_identification/` | [identify_csa_svs_by_local_overlap.R](structure_associated_SV_identification/identify_csa_svs_by_local_overlap.R) | Test associations between effective SV carrier status and local CDH occurrence. |
| `structure_associated_SV_feature/` | [analyze_csa_sv_structure_proximity.R](structure_associated_SV_feature/analyze_csa_sv_structure_proximity.R) | Compare proximity to TAD boundaries and loop anchors for CSA-SVs and other CDH-overlapping SVs. |
| `expression_analysis_mapping/` | [quantify_long_read_rna_isoquant.sh](expression_analysis_mapping/quantify_long_read_rna_isoquant.sh) | Quantify genes and transcripts from aligned long-read RNA-seq BAM files. |
| `expression_analysis_pattern/` | [test_csa_sv_expression_associations.R](expression_analysis_pattern/test_csa_sv_expression_associations.R) | Test genotype–expression associations and dosage trends. |
| `expression_analysis_integration/` | [aggregate_contact_matrices_by_sv_genotype.py](expression_analysis_integration/aggregate_contact_matrices_by_sv_genotype.py), [aggregate_insulation_by_sv_genotype.R](expression_analysis_integration/aggregate_insulation_by_sv_genotype.R) | Plot genotype-group contact matrices, pairwise matrix ratios and mean insulation profiles. |
| `expression_analysis_correlation/` | [plot_insulation_expression_correlations.R](expression_analysis_correlation/plot_insulation_expression_correlations.R) | Plot individual-level insulation–expression correlations on raw TPM and log2(TPM + 1) scales. |
| `evolution_analysis_mapping/` | [map_nhp_hic_to_chm13.sh](evolution_analysis_mapping/map_nhp_hic_to_chm13.sh), [map_nhp_hic_to_hap1.sh](evolution_analysis_mapping/map_nhp_hic_to_hap1.sh), [map_nhp_hic_to_hap2.sh](evolution_analysis_mapping/map_nhp_hic_to_hap2.sh) | Map nonhuman-primate Hi-C data, processing two replicates separately. |
| `evolution_analysis_feature/` | [compare_csa_sv_sizes_by_origin.R](evolution_analysis_feature/compare_csa_sv_sizes_by_origin.R) | Compare SV sizes using supplied ancestral/derived classifications. |

---

## Analysis workflow

1. Prepare reference assemblies, indexes, chromosome-size files and sequencing inputs; map the same human source reads to CHM13, Hap1 and Hap2.
2. Build contact matrices and call TAD boundaries and loops. Compare retained contacts and loop calls between mapping targets.
3. Prepare filtered haplotype-specific read-ID lists, extract their contact pairs, and generate matching specific and total contact counts.
4. Call CDHs separately for each sample and haplotype, then evaluate overlap with SVs and chromatin features.
5. Supply CDHs projected into CHM13 coordinates for recurrence and local SV–CDH association analyses.
6. Combine CSA-SV genotypes with prepared gene-expression tables, CHM13 contact matrices and insulation tracks for expression associations and locus-level plots.
7. Use nonhuman-primate mappings and externally prepared origin classifications for the available evolutionary analyses.

The input preparation required at each stage is documented in the corresponding script header. For example, contact extraction uses existing filtered read-ID lists, IsoQuant quantification starts from aligned BAMs, and recurrence analysis starts from projected CDH intervals.

---

## Major analysis modules

### Reference-dependent contact and loop comparisons

The mapping-comparison scripts match read-pair identifiers between CHM13 and one donor-specific haplotype. For contacts classified as cis in both mappings, the absolute contact-distance difference is `abs(D_CHM13 - D_DSA)`, in bp. The calculation uses each assembly's own positions and does not perform coordinate liftover.

To reproduce the Methods restriction to corresponding chromosomes, chromosome correspondence must be enforced before a pair contributes to the distance summary. When using the comparison logic without an internal correspondence check, supply inputs jointly filtered by read ID to satisfy this condition. Shared read IDs and separate cis filters alone do not establish chromosome correspondence. Apply this restriction to the distance-analysis inputs while preserving the input population used for any overall mapping-recovery comparison.

The distance-summary scripts count differences strictly greater than **5, 10, 25, 50 and 100 kbp**. Their `both_cis` output field counts records with a nonmissing fourth-column distance difference; this is the denominator for fractions computed from those counts. The `common_total` field counts all shared records. Wait for all comparison jobs to finish before summarizing them.

Loop comparisons are a separate analysis. `compare_loops_liftover.py` uses assembly alignments and `rustybam` to project query anchors, then classifies query loops as `BothMatch`, `SingleMatch` or `NoMatch`. Follow its header carefully for PAF direction and anchor-matching definitions.

### CDH identification and overlap analyses

`call_contact_difference_hotspots.R` models haplotype-specific contact counts using shared-contact coverage and chromosome, with a negative-binomial mean model and coverage-dependent dispersion. It applies upper-tail tests and Benjamini–Hochberg (BH) correction within each sample and haplotype. Default tested windows require total count >= 50; significant windows require q < 0.01, observed/expected >= 2 and specific count >= 10. Adjacent significant windows are merged into CDH regions.

The two SV-overlap scripts have different statistical units:

| Script | Observed count | Percentage denominator |
| --- | --- | --- |
| `test_cdh_sv_enrichment.R` | CDHs overlapping at least one SV feature interval | Number of CDHs at the selected region/window level |
| `calculate_sv_in_cdh_fraction.R` | SV records overlapping at least one CDH | Number of included SV records of the selected type |

Each counted CDH or SV contributes at most one hit, even if it overlaps multiple intervals. The first script merges SV feature intervals; the second retains individual SV records. The SV-side calculation supplies the overlap percentages summarized in the manuscript as approximately 2.9% for insertions and 4.1% for deletions. These summaries are medians of per-haplotype percentages, rather than pooled hit counts divided by pooled SV counts.

Both scripts use 1,000 permutations, repositioning CDHs on their original chromosomes while preserving interval lengths. Randomization bounds are estimated from the maximum observed CDH/SV endpoints on each chromosome. Empirical enrichment P values are reported directly. `test_cdh_structure_enrichment.R` similarly evaluates CDH overlap with TAD boundaries and loop anchors; its bounds also incorporate those feature intervals. These observed-overlap fractions and permutation-based enrichment statistics should be interpreted separately.

For haplotype-coordinate overlap analyses, Hap1 uses `hap1 INS` and `hap2 DEL`, while Hap2 uses `hap2 INS` and `hap1 DEL`. Each selected file must express the relevant sequence intervals in the coordinates of the haplotype being analyzed. This filename pairing does not perform a coordinate conversion.

### Recurrent CDHs and CSA-SV identification

`analyze_cdh_recurrence.R` uses preprojected CHM13 CDH intervals. It defines recurrent cores as continuous intervals supported by at least 20 distinct haplotypes, retaining complete cores from 12.5 kbp to 1 Mbp. The default null uses 999 circular-shift permutations, with the two haplotypes of each donor sharing the same offset on a chromosome.

`identify_csa_svs_by_local_overlap.R` tests local CDH occurrence against effective carrier status using two-sided Fisher's exact tests. The default deletion-label exchange between a donor's haplotypes is part of the effective-carrier definition. Tests use the original contingency counts; reported odds ratios add a pseudocount of 1 to each cell. BH correction is applied separately for each flank, with the default selection based on all eligible SV loci. The main `keep` set requires odds ratio > 2, adjusted P < 0.05 and at least three effective carrier haplotypes with a local CDH. The manuscript's +/-25 kbp analysis corresponds to the `f25000` outputs.

### Expression associations and chromatin integration

`test_csa_sv_expression_associations.R` uses a prepared candidate SV–gene table, with one retained observation per individual and SV–gene pair. It tests raw TPM across genotype groups using Kruskal–Wallis tests and tests dosage trends using Spearman correlation. BH correction is applied across the result table separately for the Kruskal–Wallis and dosage-trend P values (`q_kw` and `q_trend`). Candidate-gene selection and upstream expression filtering must already be reflected in the supplied table.

`aggregate_contact_matrices_by_sv_genotype.py` extracts local 25-kbp contact matrices, normalizes each individual's matrix by its total regional signal, and averages within genotype groups. With all three groups available, the output has two rows and three columns:

| Column | Top: mean normalized contact density | Bottom: log2 contact-density ratio |
| --- | --- | --- |
| 1 | 0/0 | 0/1 versus 0/0 |
| 2 | 0/1 | 1/1 versus 0/0 |
| 3 | 1/1 | 1/1 versus 0/1 |

The top row uses a shared logarithmic color scale for mean normalized density. The bottom row uses `log2((numerator_mean + eps) / (denominator_mean + eps))`, with `eps = vmax / 1000`. The script saves PDF/PNG figures and group mean matrices. Matrix values are read with `balance=False`; prepare the intended stored contact values upstream. A `.corrected.cool` filename alone does not apply Cooler balancing weights.

`aggregate_insulation_by_sv_genotype.R` reads precomputed insulation tracks and plots genotype-group mean profiles with standard-error ribbons. It does not perform pairwise genotype-group significance tests.

`plot_insulation_expression_correlations.R` produces raw-TPM and `log2(TPM + 1)` scatterplots using the same complete observations. Each point represents one individual and is colored by genotype. Plot annotations show two-sided Spearman correlations (`exact = FALSE`), nominal P values and sample sizes. Pearson results are also retained in the summary table. The regression line and confidence band describe a pooled linear fit on the displayed scale. This script applies no BH correction or covariate adjustment; its Spearman P values are distinct from the BH-adjusted genotype–expression results above.

### Evolutionary analyses

The nonhuman-primate mapping scripts generate contact-pair and alignment outputs for the supplied reference or haplotype assemblies. `compare_csa_sv_sizes_by_origin.R` takes an existing origin-classification table and compares ancestral and derived SV sizes using two-sided Wilcoxon rank-sum tests, separately for insertions and deletions. It reports nominal P values and does not infer ancestral states.

---

## Software and computational environment

Dependencies differ between modules. Read each script header and any `module load` statements before running it.

| Module | Main requirements |
| --- | --- |
| Human and nonhuman-primate mapping | Bash, BWA, SAMtools, pairtools; supporting tools listed in the scripts |
| Contact matrices and structural features | Pairix/pbgzip, cooler, HiCExplorer; Juicer Tools for the relevant matrix-conversion steps |
| Read-ID and distance comparisons | Bash, awk, gzip/zcat, GNU sort and join |
| Loop liftover | Python 3, intervaltree, rustybam |
| Haplotype-specific contacts and mappability | seqtk, BWA, pairtools, Newmap, and shell utilities as specified |
| R analyses | R with data.table, GenomicRanges, IRanges, ggplot2, MASS and splines, as required by each script |
| Matrix aggregation | Python 3 with cooler, numpy, pandas and matplotlib |
| Insulation–expression plots | R >= 4.0, data.table, ggplot2 >= 3.4.0; bash, gzip and awk for rebuild mode |
| Long-read expression quantification | IsoQuant and prepared aligned BAMs |

Versions recorded in the scripts include BWA 0.7.17, pairtools 1.0.3, SAMtools 1.19, Pairix 0.3.8, HiCExplorer 3.7.3 and IsoQuant 3.10.0. The matrix-aggregation script documents cooler 0.9.3 as the analysis version. Preserve the relevant analysis versions when reproducing results and adapt cluster-specific module names to your installation.

Many shell scripts submit Sun Grid Engine (SGE) jobs using `qsub`; others run directly. Adjust resource directives and module loading for your cluster. Submission does not imply completion of the analysis.

---

## Usage

Clone the repository:

```bash
git clone https://github.com/LingbinNi/pangenome_chromatin_associated_SVs.git
cd pangenome_chromatin_associated_SVs
```

Before running a script:

1. Read its input formats, filename patterns and coordinate conventions.
2. Replace `/path/to/...` entries and edit the settings near the beginning of the script. Most scripts use these in-file settings; scripts with command-line arguments document them explicitly.
3. Check sample names, reference assemblies, chromosome labels and required preprocessing. Use paths without whitespace where required by the shell scripts.
4. Prepare a separate output directory and follow the script's working-directory requirements. Some scripts require pre-existing reference subdirectories; others require a new output location.
5. Confirm that all upstream files and submitted jobs are complete before starting dependent analyses.

For example, after setting `cdhroot`, `svroot` and `outdir` in the SV-fraction script:

```bash
repo_dir="/path/to/pangenome_chromatin_associated_SVs"
mkdir /path/to/sv_cdh_fraction_results
cd /path/to/sv_cdh_fraction_results
Rscript "$repo_dir/contact_difference_hotspot_analysis/calculate_sv_in_cdh_fraction.R"
```

After editing the paths and settings in the matrix-aggregation script, run it from its chosen output directory:

```bash
python /path/to/pangenome_chromatin_associated_SVs/expression_analysis_integration/aggregate_contact_matrices_by_sv_genotype.py
```

The contact-counting script accepts two positional arguments:

```bash
perl /path/to/pangenome_chromatin_associated_SVs/contact_difference_hotspot_identification/count_long_range_cis_contacts.pl \
    /path/to/contact_pixels.tsv \
    /path/to/output_prefix
```

Its input is a headerless, seven-column `cooler dump --join` table containing raw integer contact counts. The parent directory for the output prefix must already exist.

Keep analysis outputs outside the code checkout. Several scripts use the current directory for output, and rerunning them can overwrite files with the same names.

---

## Reproducing manuscript analyses

The following table identifies selected calculation and plotting components. A listed script may supply statistics or panels that require separate assembly into the manuscript figure.

| Manuscript component | Script or output |
| --- | --- |
| CDH–chromatin-feature enrichment, Figure 2D | `test_cdh_structure_enrichment.R` |
| CDH-side SV-overlap enrichment, Figure 2E | `test_cdh_sv_enrichment.R`; the counted unit is a CDH |
| SV-side overlap counts, percentages and permutation summaries, Figure S5 | `calculate_sv_in_cdh_fraction.R`; the counted unit is an SV record |
| Genotype–expression association statistics, Figure 4A | `test_csa_sv_expression_associations.R`; association tables are produced, while the final expression-pattern heatmap is assembled separately |
| Genotype-group contact-density maps, Figure S9 | Top row of `aggregate_contact_matrices_by_sv_genotype.py` |
| Three genotype-group contact-map ratios, Figure 4D top | Bottom row of `aggregate_contact_matrices_by_sv_genotype.py` |
| Mean insulation profiles, Figure 4D bottom | `aggregate_insulation_by_sv_genotype.R` |
| Insulation–expression scatterplots, Figure 4E | `Insulation_Expression.log2TPM.pdf/png` from `plot_insulation_expression_correlations.R`; raw-TPM versions are also produced |
| SV-size comparisons by inferred origin, Figure S10 | `compare_csa_sv_sizes_by_origin.R`, using supplied origin labels |

For Figure 4D, the combined script output preserves the ratio order **0/1 versus 0/0, 1/1 versus 0/0, 1/1 versus 0/1**. The manuscript panel order **0/1 versus 0/0, 1/1 versus 0/1, 1/1 versus 0/0** therefore requires reordering the ratio panels during figure assembly. Insulation profiles are produced by the separate R script.

Use the sample counts, statistics and P values generated from the actual analysis inputs when preparing the manuscript. The visualization scripts do not reproduce additional significance tests merely because those tests are discussed alongside a figure.

---

## Data availability

Raw sequencing data, assemblies, contact matrices, genotype tables, insulation tracks and other cohort-scale inputs are not bundled with this repository. Input layouts and required columns are documented in the script headers. The current repository does not contain a bundled example dataset or expected-output test suite.

Data-access information and accession identifiers are described in the manuscript's **Data availability** section. Controlled-access genotype and phenotype datasets remain subject to their access requirements. Third-party software and external datasets remain subject to their respective licenses and terms of use.

---

## Reproducibility and versioning

For a reproducible run, retain the repository commit, edited script settings, input-file versions, sample lists, reference assemblies, software versions and job logs. Preserve random seeds and permutation counts where specified, and report the actual number of observations contributing to each analysis.

The `main` branch may change during manuscript preparation and revision. Record the code version with:

```bash
git rev-parse HEAD
```

The code version associated with the published manuscript will be identified in the final code-availability statement. The full cohort-scale analyses require source data and substantial computing resources.

---

## Citation

If you use code or methods from this repository, please cite:

> Ni L. *et al.* **Structural variants shape recurrent chromatin contact divergence in humans.**
> Journal information and DOI will be added upon publication.

---

## Contact

For questions regarding the code or analyses, please open a GitHub issue or contact:

**Lingbin Ni**

University of Washington, Eichler Laboratory

---

## Disclaimer

This repository contains research code developed for the analyses described in the accompanying manuscript. It is provided to support transparency and reproducibility of those analyses.
