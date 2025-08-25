---
title: "Project Outline"
date: "2025-08-24"
---

# **Section (II) Examining Two "Shot-Gun" Proteomics Experiments**
**Shotgun proteomics** is a high-throughput strategy for identifying and quantifying proteins in complex biological samples. It works by enzymatically digesting proteins into peptides, separating those peptides via liquid chromatography, and analyzing them using tandem mass spectrometry (MS/MS). This approach allows researchers to characterize thousands of proteins in a single experiment without prior knowledge of the sample’s composition.

In this project, we will explore how protein abundance information can be extracted from public datasets in the [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/projects/), using mass spectrometry-based shotgun proteomics. We focus on two biologically and technically distinct datasets:

* **Soybean (*Glycine max*)**: [PXD023343](https://www.ebi.ac.uk/pride/archive/projects/PXD023343)
* **Human (*Homo sapiens*)**: [PXD005187](https://www.ebi.ac.uk/pride/archive/projects/PXD005187)

These examples illustrate different use cases, labeling strategies, and experimental goals in proteomics.  

## Dataset Comparison: Treatment Conditions and Labeling Methods

#### **Soybean Dataset (PXD023343)**

This study investigates the response of soybean leaves to **low-phosphate (LP) stress**, a major agronomic challenge affecting crop yield. Plants were grown under **low and high phosphate conditions**, and proteins were extracted from leaves for analysis. The goal was to identify molecular pathways and proteins involved in adaptation to phosphate limitation. The experiment used **label-free quantification (LFQ)**, which estimates relative protein abundance based on peptide ion intensities without introducing chemical or metabolic labels — a practical and scalable approach for plant studies where metabolic labeling is not feasible.

#### **Human Dataset (PXD005187)**

This study explores the **molecular mechanisms by which hypoxia contributes to idiopathic pulmonary fibrosis (IPF)**, a progressive lung disease marked by epithelial injury and extracellular matrix remodeling. Human lung epithelial cells were exposed to hypoxic conditions, and shotgun proteomics was used to assess proteome-wide changes. The experiment employed **SILAC (Stable Isotope Labeling by Amino acids in Cell culture)**, a metabolic labeling technique that incorporates isotopically labeled amino acids during protein synthesis. This enables accurate MS1-level quantification of relative protein expression between normoxic and hypoxic samples. The analysis focused on hypoxia-induced shifts in signaling pathways such as **TGF-β1** and **FAK1**, and their role in promoting fibrotic phenotypes.

### Summary of Key Differences

| Feature                   | Soybean (PXD023343)                          | Human (PXD005187)                          |
| ------------------------- | -------------------------------------------- | ------------------------------------------ |
| **Biological System**     | *Glycine max* (soybean leaves)               | Human lung epithelial cells                |
| **Condition Studied**     | Low vs. high phosphate stress                | Normoxia vs. hypoxia (fibrosis model)      |
| **Quantification Method** | Label-Free Quantification (LFQ)              | SILAC (metabolic labeling)                 |
| **Study Goal**                  | Understand phosphate stress response         | Identify hypoxia-driven fibrotic signaling |

Together, these datasets demonstrate how shotgun proteomics is applied across diverse systems — from plant stress physiology to human disease — and how different labeling and quantification strategies are selected based on experimental needs and sample types.


## **(II-2): Metadata Extraction**

In this section, you'll extract important metadata from both the **original publication** and the **experimental files** associated with your dataset. This information will be used to build a `*.json*` configuration file required by SAGE to process the data.   
[1] DOI: 10.1038/cddiscovery.2017.10, Des: Kathiriya JJ, Nakra N, Nixon J, Patel PS, Vaghasiya V, Alhassani A, Tian Z, Allen-Gipson D, Davé V. Galectin-1 inhibition attenuates profibrotic signaling in hypoxia-induced pulmonary fibrosis. Cell Death Discov. 2017 Apr 10;3:17010. eCollection 2017  
[2] DOI: 10.3390/ijms22020920, Des: Cheng L, Min W, Li M, Zhou L, Hsu CC, Yang X, Jiang X, Ruan Z, Zhong Y, Wang ZY, Wang W. Quantitative Proteomics Reveals that GmENO2 Proteins Are Involved in Response to Phosphate Starvation in the Leaves of <i>Glycine max</i> L. Int J Mol Sci. 2021 22(2)  

Before editing or building your config file, read the methods section of the relevant publication and answer the following questions:

**What digestion enzyme was used and what amino acids does it cut after?**  
- *Why it matters:* The enzyme determines the peptide cleavage rules (e.g., trypsin cleaves after K or R but not before P). This affects which peptides are generated and matched during database searching, and must be correctly set in the `"enzyme"` section of `sage.conf`.
- *Where to look:* Check the methods section of the manuscript for sample preparation or digestion steps. Common entries include trypsin, Lys-C, chymotrypsin, or Glu-C. This information may also be in PRIDE metadata under “Sample Protocol” or "Digestion."

  
**Was the dataset labeled or label-free?**  
- *Why it matters:* Determines which quantification method (`tmt` vs. `lfq`) to enable in the config file.    
- *Where to look:* Look for keywords like “TMT”, “iTRAQ”, or “LFQ” in the methods section or PRIDE metadata.    
         
      
**What version of labeling was used (e.g., TMT11, TMT16, SILAC) if any?**
- *Why it matters:* Tells SAGE how many channels or mass shifts to expect (e.g., sets `"quant.tmt": "Tmt16"` or adds SILAC mods).  
- *Where to look:* Check the manuscript's methods section, figure legends, or raw file names (e.g., `TMT16`, `heavy`, `SILAC`).


**Were any post-translational modifications (PTMs) or chemical modifications used?**
- *Why it matters:* Determines what to include in `static_mods` (e.g., carbamidomethylation) and `variable_mods` (e.g., oxidation).  
- *Where to look:* Usually listed under *sample preparation* or *search parameters* in the methods section.


**Where are the `.mzML` files stored on your system?**
- *Why it matters:* You must list full or relative paths to `.mzML` files in the `"mzml_paths"` field of the config.  
- *Where to look:* The output of Step (1).

Once you have answered these questions you are ready to make a configuration file for the Open Modification Search engine known as SAGE!


## **(II-3): Open Modification Search with Spectral Alignment Guided Engine (SAGE)**
SAGE requires a configuration file (`sage.conf`) in JSON format that defines all the parameters needed for database searching, spectrum scoring, and quantification. There are many parameters where the default settings are appropriate for 90% of cases, but there are some key parameters that need to be set from the metadata you gathered above. 

This step is critical — an improperly constructed configuration file can lead to missed identifications, incorrect quantification, or complete analysis failure. An template of the configuration file is shown below. Take a moment to read through the the following sections on how SAGE works to take the raw mass spectra and match them to peptides and proteins. 
### What is SAGE?

**SAGE (Spectral Alignment Guided Engine)** is an open-source, high-throughput software tool for **database searching and quantification of mass spectrometry (MS/MS) proteomics data**. Developed by the Lazar Lab at the NIH, SAGE is designed to **scale efficiently** across large datasets while maintaining **high sensitivity and accuracy** in peptide identification.

SAGE implements modern algorithmic strategies and leverages **predictive scoring models**, **retention time prediction**, and **multi-level quantification schemes** to identify and quantify peptides from tandem mass spectra with high confidence.

SAGE is particularly useful for:

* High-throughput proteomics datasets
* Label-free and isobaric tag quantification (e.g., TMT, SILAC)
* Studies requiring open or semi-open modification searches
* Peptide identification across multiple samples with alignment

## How SAGE Works: Overview of the Pipeline
## 1. Database Generation (In Silico Digestion)

* The input FASTA file is **digested virtually** using user-defined cleavage rules (e.g., trypsin cuts after K/R but not before P).
* All possible peptide sequences are generated, filtered by:

  * Length (`min_len`, `max_len`)
  * Mass (`peptide_min_mass`, `peptide_max_mass`)
* **Static and variable modifications** are applied to produce theoretical peptide ions.
* Peptides are indexed into a search-efficient structure (controlled by `bucket_size`).

## 2. Spectral Preprocessing

* Experimental spectra (from `.mzML` files) are filtered:

  * Low-intensity noise is removed
  * Only spectra with sufficient peaks (`min_peaks`) are retained
* Optionally, deisotoping or chimeric deconvolution can be applied.

## 3. Scoring and Peptide Matching

* SAGE matches experimental MS/MS spectra to the in silico-generated peptide fragments using **scoring algorithms** based on:

  * Ion type matches (e.g., b/y ions)
  * Fragment ion intensity
  * Mass tolerance (in ppm or Da)
  * Number of matched ions (`min_matched_peaks`)
* It computes **statistical scores** to rank peptide-spectrum matches (PSMs).

## 4. Retention Time (RT) Prediction

* If enabled (`predict_rt: true`), SAGE uses machine learning to **predict the expected retention time** for each peptide.
* RT alignment helps improve scoring by matching predicted vs. observed elution behavior, reducing false positives.

## 5. Quantification

Quantification is handled differently depending on the method:

#### Label-Free Quantification (LFQ):

* SAGE integrates **MS1 precursor intensity** across chromatographic peaks.
* It scores peaks using spectral shape and alignment (`spectral_angle`, `peak_scoring`).
* Peak areas are summed or averaged (`integration: "Sum"` or `"Max"`) for protein quantification.

#### SILAC:

* Treated as a special case of LFQ, where **precursor ions are mass-shifted** due to heavy isotope incorporation.
* Peptides with the same sequence but different labels (e.g., light vs. heavy K/R) are aligned and quantified based on precursor intensity.

#### TMT (Tandem Mass Tags):

* Reporter ion intensities are extracted from **MS2 or MS3** spectra.
* Each TMT channel (e.g., TMT10, TMT16) corresponds to a different sample.
* The software optionally applies **signal-to-noise filtering** and corrects for isotopic impurities.


## 6. Output Files

* `.sage.tsv`: List of PSMs with scores
* `.quant.tsv`: Peptide-level quantification
* `.features.tsv`: Peak group features across samples
* Log files and optional diagnostic plots

### Learn More

* **GitHub**: [https://github.com/lazear/sage](https://github.com/lazear/sage)

**Defining the Clevage Enzyme Conditions**  
The enzymes used to cleage the protein samples into short peptides cut after different amino acids along the protein backbone. Some of the common digestion enzymes are listed below.  
| Enzyme           | Cleaves After (Residues)                                                     | Cleavage Blocked By | Terminal Side of Cleavage | Notes                                             |
| ---------------- | ---------------------------------------------------------------------------- | ------------------- | ------------------------- | ------------------------------------------------- |
| **Trypsin**      | K (Lysine), R (Arginine)                                                     | P (Proline)         | C-terminal                | Most commonly used; highly specific               |
| **Lys-C**        | K (Lysine)                                                                   | —                   | C-terminal                | Stable in high denaturant (e.g., urea)            |
| **Arg-C**        | R (Arginine)                                                                 | —                   | C-terminal                | Less specific than trypsin                        |
| **Glu-C**        | E (Glutamate), sometimes D (Aspartate)                                       | —                   | C-terminal                | Cleaves more broadly at high pH                   |
| **Chymotrypsin** | F (Phenylalanine), Y (Tyrosine), W (Tryptophan), L (Leucine), M (Methionine) | P (Proline)         | C-terminal                | Broader specificity; more missed cleavages likely |
| **Asp-N**        | D (Aspartate)                                                                | —                   | N-terminal                | N-terminal cleavage; less common                  |
| **Pepsin**       | Broad (especially F, L, E)                                                   | pH-dependent        | C-terminal or mixed       | Non-specific, active at low pH                    |
| **Thermolysin**  | L, I, V, F, M                                                                | —                   | N-terminal                | Metalloprotease, stable at high temperature       |


**Static vs. Variable Modifications**

In proteomics, modifications refer to changes in the mass of a residue or terminus caused by chemical labeling, post-translational modifications (PTMs), or experimental artifacts. These are defined in two categories in `sage.conf`:

**Static Modifications**

A *static modification* is a change that occurs on every instance of a given residue or terminus. The search engine always adds this mass shift to the residue — it is not optional.

Use static mods when:

* The modification was applied to all peptides during sample preparation
* You are using chemical labels like TMT or iTRAQ
* The mass shift is consistent and universal

Examples of static mods:

| Modification         | Residue | Mass Shift | When It Occurs                       |
| -------------------- | ------- | ---------- | ------------------------------------ |
| Carbamidomethylation | C       | +57.0215   | During alkylation with iodoacetamide |
| TMT16 labeling       | K, ^    | +304.207   | Chemical labeling (TMT 16plex)       |
| iTRAQ 8plex          | K, ^    | +304.2054  | iTRAQ chemical labeling              |

**Variable Modifications**

A *variable modification* is optional — the search engine will consider both the modified and unmodified versions of the peptide. This allows detection of PTMs or partial labeling.

Use variable mods when:

* You expect the modification to occur on only some peptide
* You're searching for biologically relevant PTMs
* You want to allow missed modifications (e.g., incomplete labeling)

Examples of variable mods:

| Modification                  | Residue | Mass Shift          | When It Occurs                 |
| ----------------------------- | ------- | ------------------- | ------------------------------ |
| Oxidation                     | M       | +15.9949            | Spontaneous or regulated PTM   |
| Deamidation                   | N, Q    | +0.984              | Non-enzymatic or enzymatic PTM |
| Pyro-glutamate formation      | ^Q, ^E  | −17.0265 / −18.0106 | N-terminal loss of ammonia     |
| Acetylation                   | ^       | +42.0106            | N-terminal PTM                 |
| SILAC heavy lysine (labeling) | K       | +8.0142             | Metabolic labeling             |
| SILAC heavy arginine          | R       | +10.0083            | Metabolic labeling             |

> Prefix `^` = N-terminal; `$` = C-terminal

Practical Tips

* **Do not overuse variable mods.** Each one adds **combinatorial complexity** to the search space and slows down processing.
* Use **no more than 3–5 variable mods** unless you have a strong reason.
* **Set `max_variable_mods`** in your config to limit how many mods can appear on a single peptide (e.g., 2).

Summary: When to Use Each

| Use Case                                                | Static Mod | Variable Mod |
| ------------------------------------------------------- | ---------- | ------------ |
| Carbamidomethylation on Cysteine                        | ✅ Yes      | ⛔ No         |
| TMT labeling (on K, N-term)                             | ✅ Yes      | ⛔ No         |
| Oxidation on Methionine                                 | ⛔ No       | ✅ Yes        |
| Pyro-glutamate (Q or E N-term)                          | ⛔ No       | ✅ Yes        |
| SILAC labeling (partial)                                | ⛔ No       | ✅ Yes        |
| Searching for PTMs (e.g., acetylation, phosphorylation) | ⛔ No       | ✅ Yes        |



## 2. Quantification Section
The H. Sapiens study used SILAC labeling. Other common labeling methods are also shown below.  

**TMT-Labeled Data**

TMT (Tandem Mass Tags) are chemical labels attached to peptides that allow multiplexed MS-based quantification using reporter ions. These are detected at MS2 or MS3 levels.

*Use this when:*

* The methods mention TMT10, TMT11, TMT16, or iTRAQ.
* Quantification is based on reporter ion intensities.

**SILAC-Labeled Data**

SILAC (Stable Isotope Labeling by Amino acids in Cell culture) incorporates **heavy isotopes into amino acids** (e.g., ^13C\_6^15N\_4-Arg and ^13C\_6^15N\_2-Lys) during protein synthesis, shifting precursor m/z but not fragment ions.

*Use this when:*

* The methods mention *SILAC*, *light/heavy*, or *metabolic labeling*.
* Quantification is based on *precursor intensities* for labeled peptide pairs.


> You can alternatively run **two separate searches** with different `static_mods` (one for light, one for heavy), and quantify separately.


**Label-Free Quantification (LFQ)**

LFQ uses *precursor ion intensity* without any chemical or isotopic labeling. It’s based on matching peptide features across runs.

*Use this when:*

* There is *no mention of TMT, iTRAQ, or SILAC*.
* Quantification relies on precursor ion area alone.


## **(II-6): Examine the Soybean peptide quantification file**
Let us examine the output from SAGE for the Soybean Samples. The initial info printed to screen summarizes the number of non-null elements and the datatype of each of your dataframes columns. We then display the dataframe in an easily scrollable table. 

## 1. SAGE Output File Description
1. **`peptide`**

   * **Definition:** The amino-acid sequence of the identified peptide (after in silico digestion, e.g., tryptic).
   * **Context in Sage:** Sage identifies peptides by matching MS/MS spectra to theoretical peptides from the database, including handling chimeric spectra and open modifications. ([GitHub][1], [@lazear][2])

2. **`charge`**

   * **Definition:** The precursor ion charge state observed in MS1 (e.g., +2, +3).
   * **Use:** Influences the theoretical fragmentation pattern and is used by Sage during candidate generation and scoring; different charge states of the same peptide are tracked separately in the PSM stage. ([GitHub][1])

3. **`proteins`**

   * **Definition:** The protein accession(s) the peptide maps to, in UniProt-style format (e.g., `sp|P12345|PROT_NAME`).
   * **Note:** Shared peptides (mapping to multiple proteins) introduce ambiguity; downstream aggregation logic (e.g., razor/parsimony or simple grouping) determines assignment. Sage itself reports the peptide-to-protein mapping as part of its identification output. ([GitHub][1])

4. **`q_value`**

   * **Definition:** The estimated false discovery rate (FDR)–adjusted confidence metric for the peptide-spectrum match (PSM), expressing the minimal FDR at which that identification would be accepted.
   * **Sage behavior:** Sage includes built-in FDR control (typically via target-decoy strategies and rescoring) to produce calibrated q-values so users can threshold identifications (e.g., ≤0.05) with an expected proportion of false positives. ([ResearchGate][3], [PubMed][4])

5. **`score`**

   * **Definition (in the context of Sage):** The primary peptide–spectrum match (PSM) score reflecting how well the experimental MS/MS spectrum agrees with the candidate peptide.
   * **Details:** Sage’s scoring pipeline includes initial matching and then optional **rescoring** that can incorporate additional features (e.g., retention time prediction, predicted vs observed fragment behavior, and other learned features) to improve discrimination between true and false PSMs. This composite score is used (alongside the target-decoy framework) to rank identifications before q-value estimation. ([GitHub][1], [ScienceDirect][5], [PubMed][4])
   * **Teaching note:** Because Sage supports rescoring, the `score` may reflect a learned or empirical combination of evidence beyond simple dot-product or correlation—students should inspect the particular configuration used (e.g., whether retention time prediction or other auxiliary features were enabled).

6. **`spectral_angle`**

   * **Definition:** A spectral similarity metric (often the **normalized spectral contrast angle**) between the observed experimental MS/MS spectrum and a reference/predicted spectrum.
   * **Relevance:** While not a core Sage-invented metric, spectral angle is a well-established geometric measure of similarity used in spectral matching contexts; it quantifies the angle between intensity vectors in high-dimensional fragment-ion space, with smaller angles (or higher normalized similarity) indicating closer agreement. It may be computed as part of auxiliary scoring or downstream quality assessment. ([PMC][6])
   * **Interpretation:** Because intensity vectors are normalized before angle calculation, the spectral angle is relatively robust to absolute intensity scaling; it gives a sense of how consistent the fragmentation pattern is with the expected or library-derived pattern.
   * **Teaching note:** Students can contrast this with simple dot product or correlation, and discuss why normalized geometric metrics (like spectral contrast angle) can be more sensitive to subtle spectral differences. ([PMC][6])

7. **`*.mzML`,..**

   * **Definition:** Quantitative intensity estimates (e.g., MS1 extracted ion currents or aggregated peptide-level abundances) for that peptide in the three biological/technical replicates of the **LP** condition.
   * **Context:** These intensities are the basis for relative quantification across conditions. Sage also provides quantification functionality (including LFQ-style approaches) that feeds into these, though the downstream script performs its own aggregation and normalization. ([ResearchGate][3], [@lazear][2])

## 2. View the SAGE Output File in Notebook


## **(II-7) Calculate and Examine the Protein Abundance Distributions for the Soybean Datasets**
Now that you have obtain a list of peptides that have been confidently matched to proteins of interest we need to process further to obtain estimates of the protein level abundances. \

## 1. Imports and setup
Import some useful python packages for analyzing and plotting dataframes. 

## 2. Filter the peptide tables
**Explanation:**

* Filters out low-confidence identifications: only keeps peptides with `q_value ≤ 0.05`, controlling peptide-level false discovery. This reduces noise downstream.

**Why:**
Filtering by q-value ensures you are only using peptide measurements with statistical support—propagating too many false positives would pollute protein-level inference.

**Teaching note:** The threshold (0.05) is conventional; students can try stricter (0.01) or looser as an exercise to see how discovery changes.

## 3. Select and clean relevant columns for analysis

**Explanation:**

* Constructs a subset of columns: retaining `peptide`, the original `proteins` annotation, `q_value`, and all other columns except technical metadata (`charge`, `score`, `spectral_angle`) which are not needed for quantification.
* The `proteins` field is in UniProt format like `sp|P12345|PROT_NAME;tr|P12345|PROT_NAME;...`.
* For each peptide make a row copy for every protein it could possibly map to.
* Drops any peptides that failed to parse into a valid UniProt accession.
  
**Why:**
Reduces clutter and focuses downstream computation on actual intensity data and needed identifiers.  
You need a clean, consistent protein identifier to collapse peptides into proteins. Dropping unmapped peptides avoids assigning ambiguous evidence.

**Teaching note:** Students can be asked what happens if peptides map to multiple proteins (this script doesn’t resolve shared peptides explicitly—it aggregates all evidence under each protein) and discuss how “razor” or parsimony approaches differ.

## 4. Estimate Protein Abundances
**Explanation:**

* Groups all peptides assigned to each UniProt accession and collapses their intensities per sample by taking the **median** across peptides.

**Why median?**
Median is robust to outliers (e.g., aberrant peptide measurements) and is a simple estimate of protein abundance from multiple peptide surrogates.

**Teaching note:** Alternatives include mean, weighted mean (by peptide quality), or more sophisticated models (like MaxLFQ’s ratio-based inference). Students can compare results of median vs mean as an exploration.
https://www.nature.com/articles/nmeth.3901.pdf

## 5. Compute Stats of Sample distributions

## 6. Make Box Plots of each Sample Separately

- Do you observe any noticable differences in the distribution of protein abundances between technical replicates under the same condition (i.e. LP and HP)?
- How about between conditions?
- Do these plots make sense with the statistics you calculated previously?
- What does this tell you about the proteomes of these two organisms and does this make sense for two very different organisms?

## **(II-8): Differential Protein Expression from LFQ Mass Spectrometry Data of Soybean Dataset**

*Using a robust, classical pipeline with fold-change inference via per-protein statistics*

**Explanation:**

* `welch_test`: performs Welch’s t-test comparing the LP and HP sample groups for a given protein (row of the protein-level matrix). `equal_var=False` allows unequal variances—a safer default in proteomics. `nan_policy='omit'` means missing values don’t crash the test.

**Teaching note:** Welch’s t-test is generally powerful when sample sizes are modest and distributions are approximately symmetric; the Mann–Whitney test is robust to non-normality but tests for difference in distribution location in a different way. It’s good practice to examine data (e.g., via violin/box plots) to decide if assumptions are reasonable.

## 1. Protein-level abundance statistics 

## 2. Standard Scale Renormalization of Protein Abundance

**Explanation:**

* Computes global mean and standard deviation separately for LP and HP samples (flattening all proteins × replicates).
* Z-scores each group independently: subtract mean, divide by standard deviation.

**Why:**
This rescales each condition so that its overall distribution is standardized. It removes global scale differences that might come from batch effects or loading differences.

**Caveat:**
Standardizing separately means the two groups are no longer directly on the same absolute scale—this could obscure true global shifts. Normally, for differential testing you'd either:

* Normalize all samples together (e.g., median normalization across all), or
* Use ratio-based methods where differences are preserved without groupwise rescaling.

**Teaching exercise:** Have students try skipping this step or doing a joint normalization and compare which proteins are called significant.

## 3. Statistical testing (Welch’s t-test)

**Explanation:**

* Applies Welch’s t-test per protein, comparing the standardized LP vs HP vectors (each of length 3).
* Stores the test statistic and raw p-value.

**Why Welch’s test?**
It does not assume equal variance between groups—a prudent choice for biological data where heteroskedasticity is common.

**Teaching note:** Point out that with only 3 replicates per group, degrees of freedom are low; encourage students to visualize variance and consider the effect of small sample size on power and false negatives.

## 4. Ranking and multiple testing correction

**Explanation:**

* Sorts proteins by raw significance (smallest p-value first) to prioritize candidates.
* Applies Benjamini–Hochberg FDR correction (`fdr_bh`) to control the expected proportion of false discoveries among calls.

**Why:**
Thousands of proteins may be tested; without correction, the number of false positives at p<0.05 would be unacceptably high. FDR balances discovery with error control.

**Teaching note:** Discuss differences among correction methods (`bonferroni`, `fdr_bh`, `fdr_by`, etc.) and when each is appropriate.

## 5. Extract significant proteins

**Explanation:**

* Filters proteins whose adjusted p-value is below 0.05, declaring them differentially abundant between LP and HP.
* Prints their UniProt accession IDs.

**Why:**
This yields a candidate list for biological interpretation: proteins whose standardized abundance differs beyond expected noise.

## 6. Volcano plots

A **volcano plot** is a scatter plot that simultaneously displays both **magnitude of change** and **statistical significance** for each protein (or peptide) when comparing two conditions. It’s a rapid visual filter to highlight candidates that are both substantially and confidently different.

### Axes:

* **X-axis:** Fold change (here, `Fold2AbundChange`). This is usually on a log2 scale (positive means higher in condition B, negative means higher in condition A). It reflects **effect size** — how much a protein’s abundance changes.
* **Y-axis:** $-\log_{10}(\text{adjusted p-value})$ from a statistical test (e.g., Welch’s t-test on replicate abundances). Higher values mean stronger statistical evidence against the null (i.e., more significant). Using the negative log makes small p-values (strong significance) appear at the top.

### What Each Point Represents:

Each point is a protein. Its horizontal position is how big the abundance change is, and its vertical position is how unlikely that change is due to random noise (after correcting for multiple comparisons).

### Coloring & Thresholds (as in your plot)

* **Horizontal line at adjusted p-value = 0.05:**
  Marks the significance cutoff. Points above this line are statistically significant (adj p-value < 0.05).

* **Vertical line at 0 fold change:**
  Separates proteins that increase (right) from those that decrease (left) between the two conditions.

* **Color scheme:**

  * **Bright red:** Significant *and* large change (|Fold2AbundChange| ≥ 1). These are high-confidence, high-effect candidates.
  * **Pink:** Significant but small change (|Fold2AbundChange| < 1). Statistically reliable but modest in magnitude—could be biologically meaningful if consistent across pathways.
  * **Black:** Not significant (adj p-value ≥ 0.05), regardless of fold change—these could be noise or underpowered.
