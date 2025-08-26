import os
import re
import sys
import pandas as pd
from pyteomics import fasta

################################################
def find_sage_files(directory):
    lfq_file = None
    results_file = None
    for fname in os.listdir(directory):
        print(fname)
        if fname.endswith('lfq.tsv'):
            lfq_file = os.path.join(directory, fname)
        elif fname.endswith('slim_results.sage.tsv'):
            results_file = os.path.join(directory, fname)
    return lfq_file, results_file
################################################

################################################
def combine_dfs(directory, lfq_df, results_df):

    output_file = os.path.join(directory, 'processed_sage_results.csv')
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping processing.")
        return pd.read_csv(output_file)
     
    # Merge the two DataFrames on the 'peptide' column
    charges = []
    protein_q = []
    for pep in lfq_df['peptide']:
        # print(pep)
        # Check if peptide exists in results_df
        if pep in results_df['peptide'].values:
            # print(f" - Found in results: {pep}")
            charge = results_df.loc[results_df['peptide'] == pep, 'charge'].values[0]
            # print(f"   - Charge: {charge}")
            charges.append(charge)
            protein_q.append(results_df.loc[results_df['peptide'] == pep, 'protein_q'].values[0])
        else:
            raise ValueError(f" - NOT found in results: {pep}")
    lfq_df['charge'] = charges
    lfq_df['protein_q'] = protein_q
    
    lfq_df.to_csv(output_file, index=False)
    print(f"Processed file written to: {output_file}")
    return lfq_df
################################################

################################################
def expand_by_uniprot(df):
    """
    Expands the DataFrame such that each row is duplicated for each UniProt ID
    found in the 'proteins' column, and a new column 'uniprot' is added with that ID.
    """
    rows = []
    mzML_columns = [col for col in df.columns if col.endswith('.mzML')]
    for _, row in df.iterrows():
        # Split the proteins string by ';' and extract UniProt IDs using regex
        entries = row['proteins'].split(';')
        n_entries = len(entries)
        for entry in entries:
            match = re.search(r'\|([A-Z0-9]+)\|', entry)
            if match:
                uniprot_id = match.group(1)
                new_row = row.copy()
                # for columns ending in .mzML divide by n_entries
                new_row[mzML_columns] = row[mzML_columns] / n_entries  # Divide intensity values by number of proteins
                new_row['uniprot'] = uniprot_id
                rows.append(new_row)
    
    # Combine all expanded rows into a new DataFrame
    expanded_df = pd.DataFrame(rows)
    return expanded_df
################################################

################################################
def load_fasta_sequences(fasta_path, id_regex=None):
    """Return {protein_id: sequence}. Uses first token in header unless id_regex captures a group."""
    import re
    seqs = {}
    for header, seq in fasta.read(fasta_path):
        if id_regex:
            m = re.search(id_regex, header)
            prot_id = m.group(1) if m else header.split()[0]
        else:
            prot_id = header.split()[0]
        reviewed, id, _ = prot_id.split('|')
        seqs[id] = seq
    return seqs
################################################

################################################
def protein_coverage(seqs, peptides, treat_IL_as_equal=True):
    """
    seqs: dict {protein_id: AA sequence}
    peptides: iterable of stripped peptide sequences (no mods)
    returns DataFrame with per-protein coverage stats.
    """
    def norm(s): return s.replace('I','L') if treat_IL_as_equal else s
    seqs_norm = {k: norm(v) for k, v in seqs.items()}
    peps_norm = [norm(p) for p in peptides if p and isinstance(p, str)]

    rows = []
    for pid, seq in seqs_norm.items():
        L = len(seq)
        covered = [False] * L
        for pep in peps_norm:
            start = seq.find(pep)
            while start != -1:
                for i in range(start, start + len(pep)):
                    covered[i] = True
                start = seq.find(pep, start + 1)
        cov_n = sum(covered)
        rows.append((pid, cov_n, L, cov_n / L if L else 0.0))
    return (pd.DataFrame(rows, columns=["uniprot","covered_residues","protein_length","coverage_frac"])
              .set_index("uniprot").sort_values("coverage_frac", ascending=False))
################################################
 
################################################
def main(directory, fasta_file):

    ## (1) Find SAGE files
    lfq_file, results_file = find_sage_files(directory)
    print(lfq_file, results_file)
    if not lfq_file or not results_file:
        print(f"Error: Could not find both lfq.tsv and results.sage.tsv in {directory}")
        sys.exit(1)

    print(f"Found lfq file: {lfq_file}")
    print(f"Found results file: {results_file}")

    ## (2) Load found SAGE files
    lfq_df = pd.read_csv(lfq_file, sep='\t')
    results_df = pd.read_csv(results_file, sep='\t')

    ## (3) Combine dataframes
    combined_df = combine_dfs(directory, lfq_df, results_df)
    print(combined_df)

    ## (4) Expand DataFrame by UniProt IDs
    expanded_df = expand_by_uniprot(combined_df)
    print(expanded_df)

    ## (5) strip mass offsets from peptide
    expanded_df['stripped_peptide'] = expanded_df['peptide'].str.replace(r'\[.*?\]', '', regex=True)
    print(expanded_df)

    ## (6) load the fastafile
    fasta_seqs = load_fasta_sequences(fasta_file)
    # print(fasta_seqs)

    ## (7) Compute protein coverage
    coverage_df = protein_coverage(fasta_seqs, expanded_df['stripped_peptide'].unique())
    print(coverage_df)

    ## (8) Merge coverage info back into expanded_df
    expanded_df = expanded_df.merge(coverage_df, left_on='uniprot', right_index=True, how='left')
    print(expanded_df)

    ## (9) Save the final DataFrame
    output_file = os.path.join(directory, "final_results.csv")
    expanded_df.to_csv(output_file, index=False)
    print(f"Final results saved to {output_file}")
################################################


################################################
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python CleanupSAGEresults.py <directory>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

print('NORMAL TERMINATION')