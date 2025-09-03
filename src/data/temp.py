import pandas as pd
import re
import numpy as np
# Remove modifications in parentheses
def extract_peptide(seq):
    # Removes any substring in parentheses, e.g., M(Oxidation) -> M
    return re.sub(r'\([^)]+\)', '', seq)


lfq_peptide_int = pd.read_csv('data/Soybean_OpenMS_results/lfq_peptide_intensities.tsv', sep='\t')
lfq_peptide_int['peptide'] = lfq_peptide_int['peptide_mods'].apply(extract_peptide)

print(f'LFQ Peptide Intensities:')
print(lfq_peptide_int.head())


## load the peptide to protein mapping
pep2prot_map = pd.read_csv('data/Soybean_OpenMS_results/peptide_to_protein_map.tsv')

print(f'Peptide to Protein Mapping:')
print(pep2prot_map.head())


## loop through the lfq peptides and find the proteins its assoicated with
rows = []
for idx, row in lfq_peptide_int.iterrows():
    # print(row)
    peptide = row['peptide']

    pep2prot_peptide_df = pep2prot_map[pep2prot_map['peptide'] == peptide].copy()
    associated_proteins = np.unique(pep2prot_peptide_df['protein'].tolist())
    
    pep2prot_peptide_df['proteins'] = ';'.join(associated_proteins)

    for prot in associated_proteins:
        uniprot = prot.split('|')[1]
        new_row = row.copy()
        new_row['uniprot'] = uniprot
        new_row = new_row.to_dict()
        # print(new_row, type(new_row))
        rows.append(new_row)
        # print(new_row)
    # quit()
    # uniprots = [prot.split('|')[1] for prot in associated_proteins]
    # pep2prot_peptide_df['uniprot'] = uniprots
    
    # drop the duplicate rows
    # print(f"Peptide: {peptide}, Associated Proteins: {associated_proteins}\n{pep2prot_peptide_df}")
    # quit()

    # if idx == 100:
    #     break

lfq_peptide_int = pd.DataFrame(rows)

# lfq_peptide_int['proteins'] = proteins
print(lfq_peptide_int.head())

## remove rows where proteins are empty
# lfq_peptide_int = lfq_peptide_int[lfq_peptide_int['proteins'].str.len() > 0]
# print(lfq_peptide_int.head())

# save
lfq_peptide_int.to_csv('data/Soybean_OpenMS_results/lfq_peptide_intensities_with_proteins.tsv', sep='\t', index=False)
print('Saved LFQ peptide intensities with proteins.')