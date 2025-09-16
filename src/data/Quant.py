
# --- Imports ---
import pyopenms as oms
import argparse
import os
import pandas as pd
import glob
import json
from collections import defaultdict
import csv
import subprocess

# --- OpenMSQuant class ---
class OpenMSQuant:
    def __init__(self, mzml_files, fasta_in, fasta_td, outdir, params=None):
        self.mzml_files = mzml_files
        self.fasta_in = fasta_in
        self.fasta_td = fasta_td
        self.outdir = outdir
        self.params = params or {}

    def run(self):
        # Use the original LFQ logic
        LFQ = self._lfq
        LFQ(
            self.mzml_files,
            self.fasta_in,
            self.fasta_td,
            self.outdir,
            **self.params
        )

    def _lfq(self,
        mzml_files,
        fasta_in,
        fasta_td,
        outdir,
        decoy_tag="DECOY_",
        decoy_pos="prefix",
        enzyme="Trypsin",
        missed=2,
        ppm=10.0,
        frag_da=0.02,
        fixed_mods=["Carbamidomethyl (C)"],
        var_mods=["Oxidation (M)", "Acetyl (Protein N-term)"],
        fdr_cut=0.01,
        idmap_rt_sec=5.0,
        idmap_mz_ppm=10.0):
        # ...existing code for LFQ function...
        def ensure_target_decoy_fasta(fasta, fasta_td, tag="DECOY_", position="prefix"):
            if os.path.exists(fasta_td):
                prots = []; oms.FASTAFile().load(fasta_td, prots)
                if (position == "prefix" and any(p.identifier.startswith(tag) for p in prots)) or \
                   (position == "suffix" and any(p.identifier.endswith(tag) for p in prots)):
                    return fasta_td
            targets = []; oms.FASTAFile().load(fasta, targets)
            decoy_gen = oms.DecoyGenerator()
            decoys = []
            for t in targets:
                d = oms.FASTAEntry(t)
                d.identifier = (tag + t.identifier) if position == "prefix" else (t.identifier + tag)
                aas = oms.AASequence.fromString(d.sequence)
                d.sequence = decoy_gen.reverseProtein(aas).toString()
                decoys.append(d)
            oms.FASTAFile().store(fasta_td, targets + decoys)
            return fasta_td

        def feature_finder_centroided(mzml_path):
            opts = oms.PeakFileOptions(); opts.setMSLevels([1])
            mz = oms.MSExperiment(); fh = oms.MzMLFile(); fh.setOptions(opts); fh.load(mzml_path, mz); mz.updateRanges()
            ff = oms.FeatureFinderAlgorithmPicked()
            feats = oms.FeatureMap(); seeds = oms.FeatureMap()
            params = ff.getParameters()
            ff.run(mz, feats, params, seeds)
            feats.setUniqueIds()
            feats.setPrimaryMSRunPath([mzml_path.encode()])
            return feats

        fasta_td = ensure_target_decoy_fasta(fasta_in, fasta_td, decoy_tag, decoy_pos)
        perfile_prot_ids, perfile_pep_ids = [], []
        salgo = oms.SimpleSearchEngineAlgorithm()
        p = salgo.getDefaults()
        p[b"precursor:mass_tolerance"]      = ppm
        p[b"precursor:mass_tolerance_unit"] = "ppm"
        p[b"fragment:mass_tolerance"]       = frag_da
        p[b"fragment:mass_tolerance_unit"]  = "Da"
        p[b"enzyme"]                        = enzyme
        p[b"peptide:missed_cleavages"]      = missed
        p[b"annotate:PSM"] = [b"fragment_mz_error_median_ppm", b"precursor_mz_error_ppm"]
        p[b"modifications:fixed"]    = fixed_mods
        p[b"modifications:variable"] = var_mods
        salgo.setParameters(p)

        for mzml in mzml_files:
            prot_ids, pep_ids = [], []
            salgo.search(mzml, fasta_td, prot_ids, pep_ids)
            perfile_prot_ids.append(prot_ids)
            perfile_pep_ids.append(pep_ids)

        all_prot = sum(perfile_prot_ids, [])
        all_pep  = sum(perfile_pep_ids, [])
        print(f"Identified {len(all_prot)} ProteinIdentification and {len(all_pep)} PeptideIdentification objects.")

        proteins = []
        oms.FASTAFile().load(fasta_td, proteins)

        pi = oms.PeptideIndexing()
        pi_p = pi.getDefaults()
        pi_p.setValue(b"decoy_string", decoy_tag)
        pi_p.setValue(b"decoy_string_position", decoy_pos)
        pi_p.setValue(b"enzyme:name", enzyme)
        pi_p.setValue(b"enzyme:specificity", "full")
        pi.setParameters(pi_p)

        ret = pi.run(proteins, all_prot, all_pep)
        if ret != 0:
            raise RuntimeError(f"PeptideIndexing failed with code {ret}")
        print("PeptideIndexing complete (targets + decoys present).")

        fdr = oms.FalseDiscoveryRate()
        fdr.apply(all_pep)
        print("FDR applied (q-values added).")

        for pid in all_pep:
            pid.setScoreType("q-value")
            pid.setHigherScoreBetter(False)
        print("Score type set to 'q-value' for all peptide IDs.")

        idf = oms.IDFilter()
        idf.removeEmptyIdentifications(all_pep)
        print("Removed empty identifications.")
        idf.filterHitsByScore(all_pep, float(fdr_cut))
        print("Filtered by q-value.")
        idf.removeDecoyHits(all_pep)
        print("Removed decoy identifications.")
        idf.removeEmptyIdentifications(all_pep)
        print("Removed empty identifications.")

        search_outfile = os.path.join(outdir, "searchfile_results_1perc_FDR.idXML")
        oms.IdXMLFile().store(search_outfile, all_prot, all_pep)
        print(f"Stored: {search_outfile}")

        pep2prot = defaultdict(set)
        for pid in all_pep:
            for hit in pid.getHits():
                seq = hit.getSequence().toUnmodifiedString()
                for ev in hit.getPeptideEvidences():
                    acc = ev.getProteinAccession()
                    if acc.startswith(decoy_tag):
                        continue
                    pep2prot[seq].add(acc)

        it = iter(pep2prot.items())
        for _ in range(5):
            try:
                pep, accs = next(it)
                print(pep, list(sorted(accs))[:5])
            except StopIteration:
                break

        rows = []
        for pid in all_pep:
            for hit in pid.getHits():
                seq = hit.getSequence().toUnmodifiedString()
                q   = hit.getScore()
                for ev in hit.getPeptideEvidences():
                    acc = ev.getProteinAccession()
                    if acc.startswith(decoy_tag):
                        continue
                    rows.append({
                        "peptide": seq,
                        "protein": acc,
                        "start": ev.getStart(),
                        "end": ev.getEnd(),
                        "aa_before": ev.getAABefore(),
                        "aa_after": ev.getAAAfter(),
                        "q_value": q,
                    })

        outfile = os.path.join(outdir, "peptide_to_protein_map.tsv")
        if rows:
            with open(outfile, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print(f"Wrote: {outfile}")
        else:
            print("No rows to write for peptide_to_protein_map.tsv (check filters).")

        feature_maps = [feature_finder_centroided(p) for p in mzml_files]
        run_basenames = [os.path.basename(p) for p in mzml_files]
        n_runs = len(feature_maps)
        print(f"Detected {sum(f.size() for f in feature_maps)} features across {n_runs} runs.")

        nonempty = [(i, fm) for i, fm in enumerate(feature_maps) if fm.size() > 0]
        if len(nonempty) == 0:
            raise RuntimeError("No features found in any run.")
        if len(nonempty) == 1:
            i_only, fm_only = nonempty[0]
            print(f"Only one run has features ({run_basenames[i_only]}). Using single-run LFQ.")
            spectra = oms.MSExperiment()
            oms.MzMLFile().load(mzml_files[0], spectra)
            idmap = oms.IDMapper()
            idp = idmap.getParameters()
            idp.setValue("rt_tolerance", float(idmap_rt_sec))
            idp.setValue("mz_tolerance", float(idmap_mz_ppm))
            idp.setValue("mz_measure", "ppm")
            idmap.setParameters(idp)
            idmap.annotate(fm_only, all_pep, all_prot, False, False, spectra)
            print(f"Annotated {len(all_pep)} peptide IDs onto the feature map.")
            rows = []
            for f in fm_only:
                pids = f.getPeptideIdentifications()
                if not pids:
                    continue
                hits = []
                for pid in pids:
                    hits.extend(pid.getHits())
                if not hits:
                    continue
                top = min(hits, key=lambda h: h.getScore())
                seq_mod = top.getSequence().toString()
                charge  = top.getCharge()
                intensity = float(f.getIntensity())
                rows.append([seq_mod, charge, intensity])
            out = os.path.join(outdir, "lfq_peptide_intensities.tsv")
            with open(out, "w") as fh:
                fh.write("peptide_mods\tcharge\t" + run_basenames[0] + "\n")
                for seq_mod, z, I in rows:
                    fh.write(f"{seq_mod}\t{z}\t{I:.6g}\n")
            print(f"Wrote: {out}")

        else:
            print('\nMulti-run LFQ (align → group → annotate consensus)')
            nz_idx, feature_maps = zip(*nonempty)
            feature_maps       = list(feature_maps)
            mzml_files         = [mzml_files[i] for i in nz_idx]
            run_basenames      = [run_basenames[i] for i in nz_idx]
            perfile_pep_ids    = [perfile_pep_ids[i] for i in nz_idx]
            perfile_prot_ids   = [perfile_prot_ids[i] for i in nz_idx]
            n_runs = len(feature_maps)
            idmap = oms.IDMapper()
            idp = idmap.getParameters()
            idp.setValue("rt_tolerance", float(idmap_rt_sec))
            idp.setValue("mz_tolerance", float(idmap_mz_ppm))
            idp.setValue("mz_measure", "ppm")
            idmap.setParameters(idp)
            for i, fm in enumerate(feature_maps):
                oms.IDFilter().removeEmptyIdentifications(perfile_pep_ids[i])
                exp = oms.MSExperiment()
                oms.MzMLFile().load(mzml_files[i], exp)
                idmap.annotate(fm, perfile_pep_ids[i], perfile_prot_ids[i], False, False, exp)
            print("Per-run ID mapping onto FeatureMaps complete.")
            print("Aligning feature maps…")
            ref_index = max(range(n_runs), key=lambda k: feature_maps[k].size())
            ref_fm = feature_maps[ref_index]
            aligner = oms.MapAlignmentAlgorithmPoseClustering()
            aligner.setReference(ref_fm)
            for i, fm in enumerate(feature_maps):
                if i == ref_index:
                    print(f"  Skipping alignment for reference map {i+1}/{n_runs} ({run_basenames[i]})")
                    continue
                trafo = oms.TransformationDescription()
                aligner.align(fm, trafo)
                oms.MapAlignmentTransformer().transformRetentionTimes(fm, trafo, True)
            cons_map = oms.ConsensusMap()
            grouper = oms.FeatureGroupingAlgorithmQT()
            grouper.group(feature_maps, cons_map)
            cons_map.setExperimentType("label-free")
            print(f"Linked features into {cons_map.size()} consensus features.")
            headers = cons_map.getColumnHeaders()
            for i, fm in enumerate(feature_maps):
                h = headers.get(i, oms.ColumnHeader())
                h.filename = run_basenames[i]
                h.size = fm.size()
                h.unique_id = fm.getUniqueId()
                headers[i] = h
            cons_map.setColumnHeaders(headers)
            print("Set column headers for consensus map.")
            def consensus_to_peptide_intensities(cm, n_runs):
                rows = []
                for cf in cm:
                    hits = []
                    for pid in cf.getPeptideIdentifications():
                        hits.extend(pid.getHits())
                    if not hits:
                        continue
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
            out = os.path.join(outdir, "lfq_peptide_intensities.tsv")
            with open(out, "w") as fh:
                fh.write("peptide_mods\tcharge\t" + "\t".join(run_basenames) + "\n")
                for r in rows:
                    fh.write(f"{r[0]}\t{r[1]}\t" + "\t".join(f"{x:.6g}" for x in r[2:]) + "\n")
            print(f"Wrote: {out}")

# --- SageQuant class ---
class SageQuant:
    def __init__(self, mzml_files, fasta_file, outdir, params=None):
        self.mzml_files = mzml_files
        self.fasta_file = fasta_file
        self.outdir = outdir
        self.params = params or {}

    def run(self):
        # Create SAGE control file and run quantification
        control_path = os.path.join(self.outdir, "sage_control.json")
        os.makedirs(self.outdir, exist_ok=True)
        control = self._create_control_file()
        with open(control_path, "w") as f:
            json.dump(control, f, indent=2)
        print(f"SAGE control file written to {control_path}")
     
        # Run SAGE (assume 'sage' is installed and in PATH)
        cmd = ["sage", control_path]
        print(f"Running SAGE: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"SAGE run failed: {e}")

    def _create_control_file(self):
        # Full LFQ config template for SAGE
        control = {
            "database": {
                "bucket_size": 32768,
                "enzyme": {
                    "missed_cleavages": 2,
                    "min_len": 7,
                    "max_len": 50,
                    "cleave_at": "KR",
                    "restrict": "P",
                    "c_terminal": True
                },
                "peptide_min_mass": 500.0,
                "peptide_max_mass": 5000.0,
                "ion_kinds": ["b", "y"],
                "min_ion_index": 2,
                "static_mods": {"C": 57.021464},
                "variable_mods": {"M": [15.994915]},
                "max_variable_mods": 2,
                "decoy_tag": "rev_",
                "generate_decoys": True,
                "fasta": self.fasta_file
            },
            "quant": {
                "lfq": True,
                "lfq_settings": {
                    "peak_scoring": "Hybrid",
                    "integration": "Sum",
                    "spectral_angle": 0.45,
                    "ppm_tolerance": 10.0
                }
            },
            "precursor_tol": {"ppm": [-12, 12]},
            "fragment_tol": {"ppm": [-20, 20]},
            "precursor_charge": [2, 6],
            "isotope_errors": [-1, 3],
            "deisotope": True,
            "chimera": True,
            "wide_window": False,
            "predict_rt": True,
            "min_peaks": 15,
            "max_peaks": 150,
            "min_matched_peaks": 6,
            "max_fragment_charge": 2,
            "report_psms": 1,
            "mzml_paths": self.mzml_files
        }
        # Override/add parameters from self.params
        def deep_update(d, u):
            for k, v in u.items():
                if isinstance(v, dict) and k in d and isinstance(d[k], dict):
                    deep_update(d[k], v)
                else:
                    d[k] = v
        if self.params:
            deep_update(control, self.params)
        return control

# --- main section ---
def main():
    parser = argparse.ArgumentParser(description='Quantify peptides using OpenMS or SAGE')
    parser.add_argument('--mzml_dir', required=True, help='Directory containing mzML files')
    parser.add_argument('--fasta_file', required=True, help='Path to FASTA file')
    parser.add_argument('--fasta_td_file', help='Path to target-decoy FASTA file (OpenMS only)')
    parser.add_argument('--outdir', default='Quant', help='Directory to save output files')
    parser.add_argument('--quant_type', choices=['OpenMS', 'SAGE'], default='OpenMS', help='Quantification type')
    parser.add_argument('--params', type=str, help='JSON string of additional parameters')
    args = parser.parse_args()

    mzml_files = glob.glob(os.path.join(args.mzml_dir, '*.mzML'))
    if not mzml_files:
        print(f"No mzML files found in {args.mzml_dir}")
        return
    print(f"Found {len(mzml_files)} mzML files.")

    params = json.loads(args.params) if args.params else {}

    if args.quant_type == 'OpenMS':
        if not args.fasta_td_file:
            print("fasta_td_file required for OpenMS quantification.")
            return
        quant = OpenMSQuant(mzml_files, args.fasta_file, args.fasta_td_file, args.outdir, params)
    else:
        quant = SageQuant(mzml_files, args.fasta_file, args.outdir, params)

    quant.run()

if __name__ == '__main__':
    main()
print(f'NORMAL TERMINATION')