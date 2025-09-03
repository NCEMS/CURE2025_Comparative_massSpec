import pyopenms as oms
import argparse
import os
import pandas as pd
import glob
import json
from collections import defaultdict
import csv
###############################################


####################################################

###############################################

# --- LFQ function from OpenMS_Processing.py ---
def LFQ(
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
#######################################################################################################

#######################################################################################################
def main():

    ## Get the user arguments
    parser = argparse.ArgumentParser(description='Process Peptonizer2000 results for best organism identification and Quantify the peptides')
    parser.add_argument('--input_dir', help='Directory to scan for peptonizer_result.csv files')
    parser.add_argument('--output_dir', default='Quant', help='Directory to save the output Quantification files')
    args = parser.parse_args()
    FASTAdir = 'data/FASTA/'

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)



    # Get the data/PRIDE/PXD006401/runAssessor/study_metadata.json and load it into a dictionary
    metadata_files = glob.glob(os.path.join(args.input_dir, '**', 'study_metadata.json'), recursive=True)
    print(f"Found {len(metadata_files)} study_metadata.json files.\n{metadata_files}")
    with open(metadata_files[0], 'r') as f:
        metadata = json.load(f)
    search_criteria = metadata['search_criteria']
    print(f"Loaded metadata search_criteria:\n{search_criteria}.")
    mzML_files = list(metadata['files'].keys())
    print(f"Loaded metadata mzML_files:\n{mzML_files}. {len(mzML_files)} files.")


    ## loop through the peptonizer files and get the mzML files and organisms for each 
    mzML_dict = {}
    for file in mzML_files:
        print(file)
        basename = os.path.basename(file).replace('.mzML', '')
        pep_file = [f for f in peptonizer_files if f'/{basename}_filtered70%/' in f]

        if len(pep_file) == 0:
            print(f"No peptonizer files found for {basename}.")
            continue
        if len(pep_file) > 1:
            print(f"Multiple peptonizer files found for {basename}: {pep_file}")
            quit()
        if len(pep_file) == 1:
            pep_file = pep_file[0]

        peptonizer_df = pd.read_csv(pep_file)
        # get the taxon_id for the highest score
        highest_score_row = peptonizer_df.loc[peptonizer_df['score'].astype(float).idxmax()]
        taxon_id = highest_score_row['taxon_id']
        print(f"Taxon ID for highest score in {file}: {taxon_id}")

        # Now we need to get the corresponding FASTA file
        fasta_file = os.path.join(FASTAdir, f"taxid-{taxon_id}.fasta")
        fasta_td_file = os.path.join(FASTAdir, f"taxid-{taxon_id}_td.fasta")
        if not os.path.exists(fasta_file):
            print(f"FASTA file not found for Taxon ID {taxon_id}: {fasta_file}")
            quit()
        else:
            print(f"FASTA file found for Taxon ID {taxon_id}: {fasta_file}")

        # Get the mzML file
        if taxon_id not in mzML_dict:
            mzML_dict[taxon_id] = {'mzML_files': [file], 'fasta_file': fasta_file, 'fasta_td_file': fasta_td_file}
        else:
            mzML_dict[taxon_id]['mzML_files'].append(file)

    print(mzML_dict)
    ## if no label detected
    import pyopenms as oms
    if search_criteria['labeling'] == 'none':
        print(f'{"#"*50}\nNo Labeling Detected. Using LFQ...')
        for taxon_id, files in mzML_dict.items():
            print(f"Processing {taxon_id} with files {files['mzML_files']} and FASTA {files['fasta_file']}")

            outdir = os.path.join(args.output_dir, f'taxid-{taxon_id}')
            LFQ(files['mzML_files'], files['fasta_file'], files['fasta_td_file'], outdir)
#######################################################################################################

#######################################################################################################
if __name__ == '__main__':
    main()
print(f'NORMAL TERMINATION')
#######################################################################################################