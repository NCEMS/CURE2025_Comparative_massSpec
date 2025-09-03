# OpenMS\_Processing — Label-Free Quant (LFQ) Pipeline

End-to-end **label-free quantification** in Python using **pyOpenMS**:

1. Database search with **SimpleSearchEngine** (target-decoy + PSM/peptide FDR)
2. MS1 feature detection (**FeatureFinderAlgorithmPicked**)
3. Map alignment (**PoseClustering**)
4. Feature linking to a **ConsensusMap** (QT)
5. ID→feature annotation and peptide-level intensity export
6. **Protein-level roll-up** (unique/razor/sum) — *added in this guide*

Script: `OpenMS_Processing.py`.

---

## What this pipeline does (in one picture)

```
mzML ──► SimpleSearchEngine + Target/Decoy ──► PeptideIndexing ──► FDR (q-values)
        │                                     │
        │                                     └─► peptide_to_protein_map.tsv
        └──────────────────────────────────────────────────────────────────────►

FeatureFinder (per run) ──► ID mapping ──► Alignment (PoseClustering) ──► Linking (QT)
                                                                  │
                                                                  ├─► lfq_peptide_intensities.tsv
                                                                  └─► lfq_protein_intensities.tsv  (this guide adds)
```

---

## Requirements

* Python 3.9+
* `pyopenms` (matching your OpenMS version)
* Standard libs: `argparse`, `csv`, `collections`, `os`, etc.

Tip (conda):

```bash
mamba create -n openms-lfq python=3.10 -y
mamba activate openms-lfq
pip install pyopenms
```

---

## Inputs

* **mzML files** (centroided MS1 recommended for feature picking)
* **FASTA** for the organism(s)

The script **creates or verifies** a **target+decoy** FASTA (prefix/suffix decoys via `DecoyGenerator`) before searching.&#x20;

---

## CLI usage

```bash
python OpenMS_Processing.py \
  --mzml_files sample1.mzML sample2.mzML \
  --fasta_in data/FASTA/taxid-9606.fasta \
  --fasta_td data/FASTA/taxid-9606_td.fasta \
  --enzyme Trypsin \
  --missed 2 \
  --ppm 10.0 \
  --frag_da 0.02 \
  --fixed_mods "Carbamidomethyl (C)" \
  --var_mods "Oxidation (M)" "Acetyl (Protein N-term)" \
  --fdr_cut 0.01 \
  --idmap_rt_sec 5.0 \
  --idmap_mz_ppm 10.0
```

Supported flags and defaults are defined in the `argparse` section (see lines adding `--mzml_files`, `--fasta_in`, `--fasta_td`, enzyme/tolerances/mods/FDR, etc.).&#x20;

---

## Outputs

* **`searchfile_results_1perc_FDR.idXML`** — filtered IDs (targets only; PSM/peptide q-values).&#x20;
* **`peptide_to_protein_map.tsv`** — peptide→protein positions (+ q-value) generated from peptide evidences.&#x20;
* **`lfq_peptide_intensities.tsv`** — peptide (seq+mods, charge) × runs intensity matrix:

  * Single-run: 1 column of intensities.&#x20;
  * Multi-run: N columns, one per mzML.&#x20;
* **`lfq_protein_intensities.tsv`** — **protein-level** intensities (this guide adds; see below).

---

## How it works (step-by-step)

### 1) Target+Decoy FASTA

If `--fasta_td` doesn’t already contain decoys, the script reverses each target sequence and tags it (prefix or suffix) before storing the combined FASTA.&#x20;

### 2) Search (SimpleSearchEngine)

Per-run search over target+decoy FASTA with your enzyme, missed cleavages, and mass tolerances; fixed/variable mods are configurable.&#x20;

### 3) PeptideIndexing + FDR

* Writes peptide evidences for mapping back to proteins; honors your decoy string and position.&#x20;
* Applies q-values at PSM/peptide level; sets score type to “q-value” (lower is better).&#x20;
* Removes high-q hits, decoys, and empties; stores filtered IDs to `searchfile_results_1perc_FDR.idXML`.&#x20;

### 4) Optional: peptide→protein mapping table

Exports a human-readable TSV with peptide start/end positions and neighboring residues per protein accession.&#x20;

### 5) Feature detection (per run)

**FeatureFinderAlgorithmPicked** on centroided MS1 to build FeatureMaps; script prints total feature count and handles empty maps. &#x20;

### 6) ID mapping (feature annotation)

Per-run ID→feature annotation using RT/m/z tolerances and the spectra-overload of `IDMapper::annotate`.&#x20;

### 7) Alignment & Linking

* **MapAlignmentAlgorithmPoseClustering** with an automatically chosen non-empty reference map.&#x20;
* **FeatureGroupingAlgorithmQT** builds a **ConsensusMap**; column headers carry run filenames, sizes, and IDs. &#x20;

### 8) Peptide-level LFQ export

For each consensus feature, pick the **lowest q-value** PSM, collect its sequence (with mods) and charge, then sum feature intensities per run → write `lfq_peptide_intensities.tsv`. &#x20;
(Single-run branch writes the analogous 1-column table.)&#x20;

---

## NEW: Protein-level LFQ roll-up

After writing `lfq_peptide_intensities.tsv`, we **aggregate peptides to proteins** using your in-memory **`pep2prot`** map (built from **FDR-filtered** IDs). The roll-up:

* Parses the peptide string into **unmodified** sequence for protein mapping
* Supports **shared peptide** handling policies:

  * `unique` — only peptides mapping to exactly one protein (default; cleanest)
  * `razor` — assign shared peptides to the protein with the largest current total intensity (simple MaxLFQ-like heuristic)
  * `sum` — add to **all** mapped proteins (double-counts)
* Optional `min_peptides` filter at the protein level (default 1)
* Writes **`lfq_protein_intensities.tsv`** with columns: `protein` + one column per run

You already added the function in your Helpers section and call it right after writing the peptide table in both the single-run and multi-run branches.

### Example call (multi-run branch)

```python
write_protein_lfq(
    rows,                       # peptide rows built from the ConsensusMap
    run_basenames=run_basenames,
    pep2prot=pep2prot,          # built after PeptideIndexing/FDR
    out="lfq_protein_intensities.tsv",
    shared="unique",            # or "razor"/"sum"
    min_peptides=1
)
```

### When to use which policy?

* **Start with `unique`** for conservative protein-level inferences (no shared evidence).
* **`razor`** if you want to retain some shared information without double-counting.
* **`sum`** only for exploratory totals where double-counting is acceptable.

---

## Quickstart examples

### Single run

```bash
python OpenMS_Processing.py \
  --mzml_files sample1.mzML \
  --fasta_in data/FASTA/taxid-9606.fasta
# ➜ lfq_peptide_intensities.tsv (1 column) + lfq_protein_intensities.tsv
```

(Uses the spectra-overload of `IDMapper::annotate` to attach IDs to features for that run.)&#x20;

### Multi-run

```bash
python OpenMS_Processing.py \
  --mzml_files a.mzML b.mzML c.mzML \
  --fasta_in data/FASTA/taxid-9606.fasta
# ➜ alignment (PoseClustering), linking (QT), peptide & protein LFQ across runs
```

(Reference map is the run with the most detected features; column headers are populated with run basenames.) &#x20;

---

## Tips & knobs

* **FDR cutoff**: `--fdr_cut 0.01` defaults to 1% peptide-level q-value; change if needed.&#x20;
* **Tolerances**: `--ppm`, `--frag_da`, `--idmap_rt_sec`, `--idmap_mz_ppm` control search & mapping strictness.&#x20;
* **Modifications**: Adjust `--fixed_mods` and `--var_mods` to match your protocol.&#x20;
* **Empty runs**: the script prunes entirely empty FeatureMaps before alignment/linking.&#x20;

---

## Troubleshooting

* **“No features found in any run.”**
  Ensure input mzML are **centroided** and the MS1 level is present. The script raises if all FeatureMaps are empty.&#x20;

* **ID mapping attaches few/no PSMs to features.**
  Check `--idmap_rt_sec` and `--idmap_mz_ppm`; too-tight tolerances can under-map. The per-run mapping uses the spectra-overload of `IDMapper::annotate`.&#x20;

* **Unexpected number of proteins in roll-up.**
  Try `shared="unique"` and raise `min_peptides` to 2+. If still low, inspect `peptide_to_protein_map.tsv` to verify evidences and accessions.

---

## Reproducibility & scaling

* All randomness-free; results driven by upstream tools and thresholds.
* For larger cohorts, run multiple mzML searches in parallel (e.g., GNU Parallel or job arrays) and then point the script at all mzMLs at once for alignment/linking.
* Consider writing intermediate FeatureMaps/ConsensusMaps to disk (mzML/featureXML/consensusXML) if you want to resume mid-pipeline.

---

## Extending

* **Protein inference with parsimony**: Integrate OpenMS `ProteinInference` to define minimal protein sets before roll-up.
* **PTM-aware roll-ups**: keep phosphorylated variants separate or collapse with rules by site class.
* **Intensity normalization**: add global normalization (median/quantile) prior to roll-up.

---

## File list (generated)

* `searchfile_results_1perc_FDR.idXML` — filtered identifications.&#x20;
* `peptide_to_protein_map.tsv` — peptide→protein mapping with positions and q-values.&#x20;
* `lfq_peptide_intensities.tsv` — peptide LFQ table (single or multi-run).&#x20;
* `lfq_protein_intensities.tsv` — **protein LFQ** (from roll-up in this guide).

