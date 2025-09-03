###################################################################################################################
# LFQ pipeline based on pyOpenMS user guide
# - Identification via SimpleSearch (with decoys + PSM FDR)
# - FeatureFinder (centroided)
# - Map alignment (PoseClustering; ID-free)
# - Feature linking (QT)
# - ID mapping onto consensus features
# - Exports peptide- and protein-level LFQ tables

import os
import subprocess
from collections import defaultdict
import pyopenms as oms
###################################################################################################################

###################################################################################################################
# -----------------------
# User inputs
# -----------------------

import argparse
parser = argparse.ArgumentParser(description="OpenMS LFQ pipeline")
parser.add_argument('--mzml_files', nargs='+', required=True, help='List of mzML files')
parser.add_argument('--fasta_in', required=True, help='Input FASTA file')
parser.add_argument('--fasta_td', default="data/FASTA/taxid-9606_td.fasta", help='Target+decoy FASTA file')
parser.add_argument('--decoy_tag', default="DECOY_", help='Decoy tag (prefix or suffix)')
parser.add_argument('--decoy_pos', default="prefix", choices=['prefix', 'suffix'], help='Decoy tag position')
parser.add_argument('--enzyme', default="Trypsin", help='Digestion enzyme')
parser.add_argument('--missed', type=int, default=2, help='Allowed missed cleavages')
parser.add_argument('--ppm', type=float, default=10.0, help='Precursor mass tolerance (ppm)')
parser.add_argument('--frag_da', type=float, default=0.02, help='Fragment mass tolerance (Da)')
parser.add_argument('--fixed_mods', nargs='+', default=["Carbamidomethyl (C)"], help='Fixed modifications')
parser.add_argument('--var_mods', nargs='+', default=["Oxidation (M)", "Acetyl (Protein N-term)"], help='Variable modifications')
parser.add_argument('--fdr_cut', type=float, default=0.01, help='FDR threshold')
parser.add_argument('--idmap_rt_sec', type=float, default=5.0, help='ID mapping RT tolerance (sec)')
parser.add_argument('--idmap_mz_ppm', type=float, default=10.0, help='ID mapping m/z tolerance (ppm)')
parser.add_argument('--quant', default="LFQ", help='Quantification method: LFQ or SILAC_R10K8 or SILAC_R6K4 or LYS8_ONLY')
args = parser.parse_args()

mzml_files = args.mzml_files
fasta_in   = args.fasta_in
fasta_td   = args.fasta_td
decoy_tag  = args.decoy_tag
decoy_pos  = args.decoy_pos
enzyme     = args.enzyme
missed     = args.missed
ppm        = args.ppm
frag_da    = args.frag_da
fixed_mods = args.fixed_mods
var_mods   = args.var_mods
fdr_cut    = args.fdr_cut
idmap_rt_sec = args.idmap_rt_sec
idmap_mz_ppm = args.idmap_mz_ppm
quant = args.quant
###################################################################################################################

###################################################################################################################
# ---- Presets for common metabolic labels (OpenMS/UniMod style names) ----
# Extend/modify as needed. Keys are channel names; values are sets of
# modification-name substrings that indicate the channel.
SILAC_R10K8 = {
    "light": set(),  # no SILAC mods present
    "heavy": {
        "Label:13C(6)15N(4)",  # Arg10
        "Label:13C(6)15N(2)",  # Lys8
    },
}
SILAC_R6K4 = {
    "light": set(),
    "heavy": {
        "Label:13C(6)",   # Arg6
        "Label:2H(4)",    # Lys4  (sometimes appears as 'Label:2H(4)' or with residue)
    },
}
LYS8_ONLY = {
    "light": set(),
    "heavy": {
        "Label:13C(6)15N(2)",  # Lys8
    },
}
###################################################################################################################

###################################################################################################################
# -----------------------
# Helpers
# -----------------------
def ensure_target_decoy_fasta(fasta, fasta_td, tag="DECOY_", position="prefix"):
    """Make a target+decoy FASTA using DecoyGenerator (docs recommend DECOY_ + reverse)."""
    if os.path.exists(fasta_td):
        # sanity: check if decoys exist
        prots = []; oms.FASTAFile().load(fasta_td, prots)
        if (position == "prefix" and any(p.identifier.startswith(tag) for p in prots)) or \
           (position == "suffix" and any(p.identifier.endswith(tag) for p in prots)):
            return fasta_td
    # build in-memory (per docs “More detailed example”)
    targets = []; oms.FASTAFile().load(fasta, targets)
    decoy_gen = oms.DecoyGenerator()
    decoys = []
    for t in targets:
        d = oms.FASTAEntry(t)                  # copy
        d.identifier = (tag + t.identifier) if position == "prefix" else (t.identifier + tag)
        aas = oms.AASequence.fromString(d.sequence)
        d.sequence = decoy_gen.reverseProtein(aas).toString()   # reverse
        decoys.append(d)
    oms.FASTAFile().store(fasta_td, targets + decoys)
    return fasta_td

def feature_finder_centroided(mzml_path):
    """FeatureFinderAlgorithmPicked on centroided data (guide shows MS1-only loading)."""
    opts = oms.PeakFileOptions(); opts.setMSLevels([1])
    mz = oms.MSExperiment(); fh = oms.MzMLFile(); fh.setOptions(opts); fh.load(mzml_path, mz); mz.updateRanges()
    ff = oms.FeatureFinderAlgorithmPicked()
    feats = oms.FeatureMap(); seeds = oms.FeatureMap()
    params = ff.getParameters()   # tune if needed
    ff.run(mz, feats, params, seeds)
    feats.setUniqueIds()
    # keep a breadcrumb to input file (useful when filling column headers)
    feats.setPrimaryMSRunPath([mzml_path.encode()])
    return feats

def view_peptide_ids(peptide_ids, fdr_cut):
    sig_count = 0
    for peptide_id in peptide_ids:
        # Peptide identification values
        print(35 * "=")
        print("Peptide ID m/z:", peptide_id.getMZ())
        print("Peptide ID rt:", peptide_id.getRT())
        print("Peptide scan index:", peptide_id.getMetaValue("scan_index"))
        print("Peptide scan name:", peptide_id.getMetaValue("scan_index"))
        print("Peptide ID score type:", peptide_id.getScoreType())
        num_hits = len(peptide_id.getHits())
        print("Peptide ID num hits:", num_hits)
        # PeptideHits
        for hit in peptide_id.getHits():
            print(" - Peptide hit rank:", hit.getRank())
            print(" - Peptide hit charge:", hit.getCharge())
            print(" - Peptide hit sequence:", hit.getSequence())
            mz = (
                hit.getSequence().getMonoWeight(
                    oms.Residue.ResidueType.Full, hit.getCharge()
                )
                / hit.getCharge()
            )
            print(" - Peptide hit monoisotopic m/z:", mz)
            print(
                " - Peptide ppm error:", abs(mz - peptide_id.getMZ()) / mz * 10**6
            )
            print(" - Peptide hit score:", hit.getScore())
            if hit.getScore() <= fdr_cut:
                sig_count += 1
    print(f"Significant peptide IDs (FDR < {fdr_cut}): {sig_count}")

def write_protein_lfq(
    rows,
    run_basenames,
    pep2prot,
    out="lfq_protein_intensities.tsv",
    shared="unique",        # "unique" (default), "razor", or "sum"
    min_peptides=1          # minimum #peptides contributing to keep a protein
):
    """
    Roll up peptide-level LFQ rows to protein-level intensities.

    Parameters
    ----------
    rows : list
        Single-run rows:    [seq_mod, charge, I]
        Multi-run rows:     [seq_mod, charge, I1, I2, ...]
    run_basenames : list[str]
        Column headers (one per run).
    pep2prot : dict[str, set[str]]
        Mapping: unmodified peptide sequence -> set of protein accessions.
        (Built earlier from FDR-filtered identifications.)
    out : str
        Output TSV path for protein-level intensities.
    shared : str
        How to handle peptides mapping to multiple proteins:
          - "unique": use only peptides that map to exactly one protein
          - "razor": assign shared peptide to the protein currently having the
                     largest accumulated intensity (simple heuristic)
          - "sum":   add the peptide to all mapped proteins (double-counts)
    min_peptides : int
        Require at least this many contributing peptides per protein.
    """
    n_runs = len(run_basenames)
    prot_int = defaultdict(lambda: [0.0] * n_runs)
    prot_pep_counts = defaultdict(int)

    def intens_vector(r):
        # r = [seq_mod, charge, ...intensities...]
        return [float(r[2])] if len(r) == 3 else [float(x) for x in r[2:]]

    for r in rows:
        seq_mod = r[0]
        intens  = intens_vector(r)

        # Convert modified string to unmodified for mapping
        try:
            unmod = oms.AASequence.fromString(seq_mod).toUnmodifiedString()
        except Exception:
            # If parsing fails, skip this peptide
            continue

        prots = list(pep2prot.get(unmod, []))
        if not prots:
            continue

        if shared == "unique":
            if len(prots) != 1:
                continue
            targets = prots
        elif shared == "sum":
            targets = prots
        elif shared == "razor":
            # assign to the protein with the largest current total intensity
            def current_total(p):
                return sum(prot_int[p])
            best = max(prots, key=current_total)
            targets = [best]
        else:
            raise ValueError("shared must be one of {'unique','razor','sum'}")

        for p in targets:
            prot_pep_counts[p] += 1
            acc = prot_int[p]
            for i in range(n_runs):
                if i < len(intens):
                    acc[i] += intens[i]

    # Filter proteins with insufficient peptide evidence
    kept = {p: vec for p, vec in prot_int.items() if prot_pep_counts[p] >= min_peptides}

    with open(out, "w") as fh:
        fh.write("protein\t" + "\t".join(run_basenames) + "\n")
        for p, vec in kept.items():
            fh.write(p + "\t" + "\t".join(f"{x:.6g}" for x in vec) + "\n")

    print(f"Wrote: {out} ({len(kept)} proteins; shared={shared}; min_peptides={min_peptides})")

def _mods_in_sequence(seq_mod_str):
    """Return a set of modification names present in AASequence (residue-level)."""
    aas = oms.AASequence.fromString(seq_mod_str)
    mods = set()
    for i in range(aas.size()):
        if aas[i].isModified():
            try:
                mods.add(aas[i].getModificationName())
            except Exception:
                # Some builds require retrieving via residue->getModification()
                m = aas[i].getModification()
                if m:
                    mods.add(m.getFullId())
    # N-term / C-term labels (rare in SILAC but common in dimethyl) 
    if aas.hasNTerminalModification():
        mods.add(aas.getNTerminalModificationName())
    if aas.hasCTerminalModification():
        mods.add(aas.getCTerminalModificationName())
    return mods

def _channel_of(seq_mod_str, channel_map):
    """
    Decide which channel a peptide belongs to based on its modification set.
    - channel_map: dict {channel -> set of substrings}; empty set means "no labels"
    Returns channel or None if ambiguous/unmatched.
    """
    mods = _mods_in_sequence(seq_mod_str)
    # Normalize to string list and match by substring (robust to residue suffixes)
    mod_texts = list(mods)
    matches = []
    for ch, needles in channel_map.items():
        if not needles:  # 'light' channel
            # belongs to 'light' only if no channel-defining labels present at all
            # i.e., if no other channel matches
            continue
        if any(any(needle in m for m in mod_texts) for needle in needles):
            matches.append(ch)
    if matches:
        # If multiple channels match (shouldn't for SILAC), declare ambiguous
        return matches[0] if len(matches) == 1 else None
    # If no heavy/mid channels matched, check for 'light' definition
    return "light" if "light" in channel_map else None

def quantify_metabolic_labels(
    rows,
    run_basenames,
    pep2prot,
    channel_map=SILAC_R10K8,
    out_prefix="silac",
    ratio_pairs=(("heavy", "light"),),   # list of (numerator, denominator)
    min_peptides_protein=1,
    protein_shared="unique",             # reuse your policy: 'unique' | 'razor' | 'sum'
):
    """
    Build channel-wise intensities and ratios for metabolic labels (e.g., SILAC).

    Parameters
    ----------
    rows : list
        From your LFQ branch:
         - single-run: [seq_mod, charge, I]
         - multi-run:  [seq_mod, charge, I1, I2, ...]
    run_basenames : list[str]
        Per-run column names.
    pep2prot : dict[str, set[str]]
        Unmodified peptide -> protein accessions (FDR-filtered) mapping.
    channel_map : dict[str, set[str]]
        e.g., {'light': set(), 'heavy': {'Label:13C(6)15N(4)', 'Label:13C(6)15N(2)'}}
    out_prefix : str
        Prefix for output files.
    ratio_pairs : list[tuple[str, str]]
        Ratios to compute, e.g., [('heavy','light')].
    min_peptides_protein : int
        Min contributing peptides to keep a protein in protein tables.
    protein_shared : str
        How to handle shared peptides in protein roll-up.
    """
    n_runs = len(run_basenames)
    # peptide key = (unmodified sequence, charge)
    peptab = defaultdict(lambda: {ch: [0.0] * n_runs for ch in channel_map.keys()})

    def intens_vector(r):
        return [float(r[2])] if len(r) == 3 else [float(x) for x in r[2:]]

    # --- Fill channel buckets at peptide level ---
    for r in rows:
        seq_mod = r[0]
        z       = int(r[1])
        I       = intens_vector(r)

        channel = _channel_of(seq_mod, channel_map)
        if channel is None:
            # ambiguous or unlabeled in a scheme that doesn't define 'light'
            continue

        try:
            unmod = oms.AASequence.fromString(seq_mod).toUnmodifiedString()
        except Exception:
            continue

        acc = peptab[(unmod, z)][channel]
        for i in range(min(n_runs, len(I))):
            acc[i] += I[i]

    # --- Write channel-wise peptide intensities ---
    pep_int_file = f"{out_prefix}_peptide_channel_intensities.tsv"
    with open(pep_int_file, "w") as fh:
        header = ["peptide_unmod", "charge"]
        for ch in channel_map.keys():
            header.extend([f"{ch}:{rb}" for rb in run_basenames])
        fh.write("\t".join(header) + "\n")
        for (pep, z), ch2vec in peptab.items():
            row = [pep, str(z)]
            for ch in channel_map.keys():
                row.extend(f"{x:.6g}" for x in ch2vec[ch])
            fh.write("\t".join(row) + "\n")
    print(f"Wrote: {pep_int_file}")

    # --- Write peptide-level channel ratios ---
    pep_ratio_file = f"{out_prefix}_peptide_ratios.tsv"
    with open(pep_ratio_file, "w") as fh:
        header = ["peptide_unmod", "charge"]
        for num, den in ratio_pairs:
            header.extend([f"{num}/{den}:{rb}" for rb in run_basenames])
        fh.write("\t".join(header) + "\n")
        for (pep, z), ch2vec in peptab.items():
            row = [pep, str(z)]
            for num, den in ratio_pairs:
                num_v = ch2vec.get(num, [0.0] * n_runs)
                den_v = ch2vec.get(den, [0.0] * n_runs)
                for i in range(n_runs):
                    nv, dv = float(num_v[i]), float(den_v[i])
                    ratio = (nv / dv) if dv > 0 else float("nan")
                    row.append(f"{ratio:.6g}")
            fh.write("\t".join(row) + "\n")
    print(f"Wrote: {pep_ratio_file}")

    # --- Protein-level per-channel roll-up (reuse your write_protein_lfq) ---
    # Build rows_like for each channel, so we can call write_protein_lfq() separately per channel.
    channel_to_rows = {ch: [] for ch in channel_map.keys()}
    for (pep, z), ch2vec in peptab.items():
        for ch, vec in ch2vec.items():
            # Reconstruct a "row": [seq_mod, charge, I...] where seq_mod is the UNMOD seq
            # (Protein roll-up ignores mods anyway via toUnmodifiedString)
            channel_to_rows[ch].append([pep, z] + vec)

    prot_tables = {}
    for ch, ch_rows in channel_to_rows.items():
        out = f"{out_prefix}_protein_{ch}_intensities.tsv"
        # NOTE: write_protein_lfq expects rows like [seq_mod, charge, I...] and will
        # parse AASequence; since we pass unmodified strings, parsing is still valid.
        write_protein_lfq(
            rows=ch_rows,
            run_basenames=run_basenames,
            pep2prot=pep2prot,
            out=out,
            shared=protein_shared,
            min_peptides=min_peptides_protein
        )
        prot_tables[ch] = out

    # --- Protein-level channel ratios ---
    # Simple merge of the per-channel protein tables on protein accession.
    # (No external deps; do an in-memory join.)
    # Load per-channel tables:
    ch_to_prot = {}
    proteins_all = set()
    for ch, path in prot_tables.items():
        table = {}
        with open(path) as fh:
            header = next(fh).rstrip("\n").split("\t")
            cols = header[1:]  # run columns
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                prot  = parts[0]
                vals  = [float(x) if x and x.lower() != "nan" else float("nan") for x in parts[1:]]
                table[prot] = vals
                proteins_all.add(prot)
        ch_to_prot[ch] = (cols, table)

    prot_ratio_file = f"{out_prefix}_protein_ratios.tsv"
    with open(prot_ratio_file, "w") as fh:
        header = ["protein"]
        for num, den in ratio_pairs:
            header.extend([f"{num}/{den}:{rb}" for rb in run_basenames])
        fh.write("\t".join(header) + "\n")
        for prot in sorted(proteins_all):
            row = [prot]
            for num, den in ratio_pairs:
                _, num_tab = ch_to_prot.get(num, (run_basenames, {}))
                _, den_tab = ch_to_prot.get(den, (run_basenames, {}))
                num_v = num_tab.get(prot, [float("nan")] * n_runs)
                den_v = den_tab.get(prot, [float("nan")] * n_runs)
                for i in range(n_runs):
                    nv, dv = num_v[i], den_v[i]
                    ratio = (nv / dv) if (dv and dv > 0) else float("nan")
                    row.append(f"{ratio:.6g}")
            fh.write("\t".join(row) + "\n")
    print(f"Wrote: {prot_ratio_file}")
###################################################################################################################

###################################################################################################################
# -----------------------
# Identification (SimpleSearch)  — docs: basic example + decoys/FDR
# -----------------------
# 1) ensure target+decoy DB
fasta_td = ensure_target_decoy_fasta(fasta_in, fasta_td, decoy_tag, decoy_pos)

# 2) search each run
perfile_prot_ids, perfile_pep_ids = [], []
salgo = oms.SimpleSearchEngineAlgorithm()
p = salgo.getDefaults()
p[b"precursor:mass_tolerance"]      = ppm
p[b"precursor:mass_tolerance_unit"] = "ppm"
p[b"fragment:mass_tolerance"]       = frag_da
p[b"fragment:mass_tolerance_unit"]  = "Da"
p[b"enzyme"]                        = enzyme
p[b"peptide:missed_cleavages"]      = missed
# Optional PSM annotations (from guide)
p[b"annotate:PSM"] = [b"fragment_mz_error_median_ppm", b"precursor_mz_error_ppm"]
p[b"modifications:fixed"]    = fixed_mods
p[b"modifications:variable"] = var_mods
salgo.setParameters(p)

for mzml in mzml_files:
    prot_ids, pep_ids = [], []
    salgo.search(mzml, fasta_td, prot_ids, pep_ids)
    perfile_prot_ids.append(prot_ids)
    perfile_pep_ids.append(pep_ids)

# 3) merge IDs across runs
all_prot = sum(perfile_prot_ids, [])
all_pep  = sum(perfile_pep_ids, [])
print(f"Identified {len(all_prot)} ProteinIdentification and {len(all_pep)} PeptideIdentification objects.")

# --- Peptide↔Protein indexing (creates peptide evidences used for mapping) ---
proteins = []
oms.FASTAFile().load(fasta_td, proteins)

pi = oms.PeptideIndexing()
pi_p = pi.getDefaults()
pi_p.setValue(b"decoy_string", decoy_tag)              # e.g., "DECOY_"
pi_p.setValue(b"decoy_string_position", decoy_pos)     # "prefix" or "suffix"
pi_p.setValue(b"enzyme:name", enzyme)
pi_p.setValue(b"enzyme:specificity", "full")
pi.setParameters(pi_p)

ret = pi.run(proteins, all_prot, all_pep)
if ret != 0:
    raise RuntimeError(f"PeptideIndexing failed with code {ret}")
print("PeptideIndexing complete (targets + decoys present).")

# --- FDR on PSM/peptide IDs ---
fdr = oms.FalseDiscoveryRate()
fdr.apply(all_pep)  # adds q-values at PSM/peptide level
print("FDR applied (q-values added).")

# make score intent explicit (helps some builds)
for pid in all_pep:
    pid.setScoreType("q-value")
    pid.setHigherScoreBetter(False)
print("Score type set to 'q-value' for all peptide IDs.")

# --- Filter by q-value; drop decoys & empty ---
idf = oms.IDFilter()
idf.removeEmptyIdentifications(all_pep)               # pre-clean to avoid crashes on some builds
print("Removed empty identifications.")
idf.filterHitsByScore(all_pep, float(fdr_cut))  # keep q ≤ fdr_cut  (True = lower-better)
print("Filtered by q-value.")
idf.removeDecoyHits(all_pep)
print("Removed decoy identifications.")
idf.removeEmptyIdentifications(all_pep)
print("Removed empty identifications.")

# Store filtered IDs (optional)
oms.IdXMLFile().store("searchfile_results_1perc_FDR.idXML", all_prot, all_pep)
print("Stored: searchfile_results_1perc_FDR.idXML")

# --- Build peptide → proteins mapping (uses peptide evidences written by PeptideIndexing) ---
from collections import defaultdict
pep2prot = defaultdict(set)

for pid in all_pep:
    for hit in pid.getHits():
        seq = hit.getSequence().toUnmodifiedString()
        for ev in hit.getPeptideEvidences():
            acc = ev.getProteinAccession()
            if acc.startswith(decoy_tag):  # skip decoys if desired
                continue
            pep2prot[seq].add(acc)

# Example: peek a few mappings
it = iter(pep2prot.items())
for _ in range(5):
    try:
        pep, accs = next(it)
        print(pep, list(sorted(accs))[:5])
    except StopIteration:
        break

# --- (Optional) write detailed peptide→protein positions with q-values ---
import csv
rows = []
for pid in all_pep:
    for hit in pid.getHits():
        seq = hit.getSequence().toUnmodifiedString()
        q   = hit.getScore()  # q-value post-FDR
        for ev in hit.getPeptideEvidences():
            acc = ev.getProteinAccession()
            if acc.startswith(decoy_tag):
                continue
            rows.append({
                "peptide": seq,
                "protein": acc,
                "start": ev.getStart(),         # 1-based positions
                "end": ev.getEnd(),
                "aa_before": ev.getAABefore(),
                "aa_after": ev.getAAAfter(),
                "q_value": q,
            })

if rows:
    with open("peptide_to_protein_map.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("Wrote: peptide_to_protein_map.tsv")
else:
    print("No rows to write for peptide_to_protein_map.tsv (check filters).")
###################################################################################################################

###################################################################################################################
# ============================================================
# LFQ intensities for the filtered peptides (per run)
# ============================================================
# 1) Feature detection per run
# for mzml_f in mzml_files:
#     fm = feature_finder_centroided(mzml_f)

feature_maps = [feature_finder_centroided(p) for p in mzml_files]
run_basenames = [os.path.basename(p) for p in mzml_files]
n_runs = len(feature_maps)
print(f"Detected {sum(f.size() for f in feature_maps)} features across {n_runs} runs.")

# Remove truly empty maps to avoid crashes
nonempty = [(i, fm) for i, fm in enumerate(feature_maps) if fm.size() > 0]
if len(nonempty) == 0:
    raise RuntimeError("No features found in any run.")
if len(nonempty) == 1:
    # fall back to your single-run branch (no alignment/linking)
    i_only, fm_only = nonempty[0]
    print(f"Only one run has features ({run_basenames[i_only]}). Using single-run LFQ.")
    # ... call your single-run annotate→export code here ...
else:
    # Rebuild arrays with only non-empty maps (keeps indexes compact 0..N-1)
    print(f"Keeping {len(nonempty)} non-empty feature maps for multi-run LFQ.")

# Common ID mapping parameters
idmap = oms.IDMapper()
idp = idmap.getParameters()
idp.setValue("rt_tolerance", float(idmap_rt_sec))   # seconds
idp.setValue("mz_tolerance", float(idmap_mz_ppm))   # ppm
idp.setValue("mz_measure", "ppm")
idmap.setParameters(idp)

if n_runs == 1:
    print('Single-run LFQ (no grouping needed)')
    # -------- Single-run LFQ (no grouping needed) --------
    fm = feature_maps[0]
    print(f"FeatureMap size: {fm.size()}")

    # --- load the raw spectra for this mzML (needed for the annotate(.., .., .., .., .., spectra) overload) ---
    spectra = oms.MSExperiment()
    oms.MzMLFile().load(mzml_files[0], spectra)

    # --- IMPORTANT: call the 6-argument overload as in the docs ---
    # clear_ids=False (keep any existing), map_ms1=False (don't try to map MS1 scans)
    idmap.annotate(fm, all_pep, all_prot, False, False, spectra)
    print(f"Annotated {len(all_pep)} peptide IDs onto the feature map.")

    # Build peptide (sequence+mods, charge) intensity table from FeatureMap
    rows = []
    for f in fm:
        pids = f.getPeptideIdentifications()
        if not pids:
            continue
        hits = []
        for pid in pids:
            hits.extend(pid.getHits())
        if not hits:
            continue
        # choose lowest q-value (we set q-values & lower-better above)
        top = min(hits, key=lambda h: h.getScore())
        seq_mod = top.getSequence().toString()   # includes modifications
        charge  = top.getCharge()
        intensity = float(f.getIntensity())      # MS1 feature intensity
        rows.append([seq_mod, charge, intensity])

    # Write final table: one column for this mzML
    out = "lfq_peptide_intensities.tsv"
    with open(out, "w") as fh:
        fh.write("peptide_mods\tcharge\t" + run_basenames[0] + "\n")
        for seq_mod, z, I in rows:
            fh.write(f"{seq_mod}\t{z}\t{I:.6g}\n")
    print(f"Wrote: {out}")


    # --- SILAC / metabolic label quantification (single run) ---
    if quant != "LFQ":
        quantify_metabolic_labels(
            rows=rows,
            run_basenames=[run_basenames[0]],
            pep2prot=pep2prot,
            channel_map=SILAC_R10K8,
            out_prefix="silac",
            ratio_pairs=(("heavy", "light"),),
            min_peptides_protein=1,
            protein_shared="unique"
        )

    # ============================================================
    # Protein intensities 
    # ============================================================
    # After: print(f"Wrote: {out}")
    write_protein_lfq(
        rows,
        run_basenames=[run_basenames[0]],
        pep2prot=pep2prot,
        out="lfq_protein_intensities.tsv",
        shared="unique",       # or "razor"/"sum"
        min_peptides=1
    )

else:
    print('Multi-run LFQ (align → group → annotate consensus)')

    # -------- keep only non-empty feature maps and sync all arrays --------
    nz_idx, feature_maps = zip(*nonempty)
    feature_maps       = list(feature_maps)
    mzml_files         = [mzml_files[i] for i in nz_idx]
    run_basenames      = [run_basenames[i] for i in nz_idx]
    perfile_pep_ids    = [perfile_pep_ids[i] for i in nz_idx]
    perfile_prot_ids   = [perfile_prot_ids[i] for i in nz_idx]
    n_runs = len(feature_maps)

    # -------- ID mapping per run (FeatureMap overload with spectra) --------
    idmap = oms.IDMapper()
    idp = idmap.getParameters()
    idp.setValue("rt_tolerance", float(idmap_rt_sec))   # seconds
    idp.setValue("mz_tolerance", float(idmap_mz_ppm))   # ppm
    idp.setValue("mz_measure", "ppm")
    idmap.setParameters(idp)

    for i, fm in enumerate(feature_maps):
        # guard: drop empty PSM containers to avoid crashes
        oms.IDFilter().removeEmptyIdentifications(perfile_pep_ids[i])

        # load spectra for THIS run only
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_files[i], exp)

        # annotate features with (already FDR-filtered) IDs from the same run
        # clear_ids=False (keep any existing), map_ms1=False
        idmap.annotate(fm, perfile_pep_ids[i], perfile_prot_ids[i], False, False, exp)

    print("Per-run ID mapping onto FeatureMaps complete.")

    # -------- Alignment (pose clustering) --------
    print("Aligning feature maps…")
    ref_index = max(range(n_runs), key=lambda k: feature_maps[k].size())
    ref_fm = feature_maps[ref_index]

    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    aligner.setReference(ref_fm)  # IMPORTANT: set a non-empty reference

    for i, fm in enumerate(feature_maps):
        if i == ref_index:
            print(f"  Skipping alignment for reference map {i+1}/{n_runs} ({run_basenames[i]})")
            continue
        trafo = oms.TransformationDescription()
        aligner.align(fm, trafo)
        oms.MapAlignmentTransformer().transformRetentionTimes(fm, trafo, True)

    # -------- Linking → ConsensusMap --------
    cons_map = oms.ConsensusMap()
    grouper = oms.FeatureGroupingAlgorithmQT()
    grouper.group(feature_maps, cons_map)
    cons_map.setExperimentType("label-free")
    print(f"Linked features into {cons_map.size()} consensus features.")

    # fill column headers
    headers = cons_map.getColumnHeaders()
    for i, fm in enumerate(feature_maps):
        h = headers.get(i, oms.ColumnHeader())
        h.filename = run_basenames[i]
        h.size = fm.size()
        h.unique_id = fm.getUniqueId()
        headers[i] = h
    cons_map.setColumnHeaders(headers)
    print("Set column headers for consensus map.")

    # -------- (Optional) extra mapping on ConsensusMap --------
    # Not strictly needed because IDs attached to FeatureMaps are propagated to the ConsensusMap.
    # If you still want to run it, use the 3-arg overload (no spectra):
    # idmap.annotate(cons_map, all_pep, all_prot)

    # -------- Build peptide (seq+mods, charge) × runs intensity table --------
    def consensus_to_peptide_intensities(cm, n_runs):
        rows = []
        for cf in cm:
            hits = []
            for pid in cf.getPeptideIdentifications():
                hits.extend(pid.getHits())
            if not hits:
                continue
            # we set q-values (lower is better), so pick the lowest q
            top = min(hits, key=lambda h: h.getScore())
            seq_mod = top.getSequence().toString()
            charge  = top.getCharge()

            idx_to_int = {}
            for h in cf.getFeatureList():
                idx = h.getMapIndex()
                idx_to_int[idx] = idx_to_int.get(idx, 0.0) + float(h.getIntensity())

            intensities = [idx_to_int.get(i, 0.0) for i in range(n_runs)]
            rows.append([seq_mod, charge] + intensities)
        return rows

    rows = consensus_to_peptide_intensities(cons_map, n_runs)

    with open("lfq_peptide_intensities.tsv", "w") as fh:
        fh.write("peptide_mods\tcharge\t" + "\t".join(run_basenames) + "\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t" + "\t".join(f"{x:.6g}" for x in r[2:]) + "\n")
    print("Wrote: lfq_peptide_intensities.tsv")

    # --- SILAC / metabolic label quantification ---
    if quant != "LFQ":
        quantify_metabolic_labels(
            rows=rows,
            run_basenames=run_basenames,
            pep2prot=pep2prot,
            channel_map=SILAC_R10K8,           # or SILAC_R6K4, LYS8_ONLY, or your custom map
            out_prefix="silac",
            ratio_pairs=(("heavy", "light"),), # defines which ratios to compute
            min_peptides_protein=1,
            protein_shared="unique"
        )

    # ============================================================
    # Protein intensities 
    # ============================================================
    # After: print(f"Wrote: {out}")
    write_protein_lfq(
        rows,
        run_basenames=run_basenames,
        pep2prot=pep2prot,
        out="lfq_protein_intensities.tsv",
        shared="unique",       # or "razor"/"sum"
        min_peptides=1
    )
###################################################################################################################

