import os
import sys
import pandas as pd

def find_sage_files(directory):
    lfq_file = None
    results_file = None
    for fname in os.listdir(directory):
        if fname.endswith('lfq.tsv'):
            lfq_file = os.path.join(directory, fname)
        elif fname.endswith('results.sage.tsv'):
            results_file = os.path.join(directory, fname)
    return lfq_file, results_file

def main(directory):
    lfq_file, results_file = find_sage_files(directory)
    if not lfq_file or not results_file:
        print(f"Error: Could not find both lfq.tsv and results.sage.tsv in {directory}")
        sys.exit(1)

    print(f"Found lfq file: {lfq_file}")
    print(f"Found results file: {results_file}")

    lfq_df = pd.read_csv(lfq_file, sep='\t')
    results_df = pd.read_csv(results_file, sep='\t')

    print(lfq_df)
    charges = []
    for pep in lfq_df['peptide']:
        # print(pep)
        # Check if peptide exists in results_df
        if pep in results_df['peptide'].values:
            # print(f" - Found in results: {pep}")
            charge = results_df.loc[results_df['peptide'] == pep, 'charge'].values[0]
            # print(f"   - Charge: {charge}")
            charges.append(charge)
        else:
            raise ValueError(f" - NOT found in results: {pep}")
    lfq_df['charge'] = charges
    print(lfq_df)

    output_file = os.path.join(directory, 'processed_sage_results.tsv')
    lfq_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed file written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python CleanupSAGEresults.py <directory>")
        sys.exit(1)
    main(sys.argv[1])

print('NORMAL TERMINATION')