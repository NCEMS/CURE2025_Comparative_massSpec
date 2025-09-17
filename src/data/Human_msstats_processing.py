import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
from pyteomics import fasta

############################################################################
def coverage(df, threshold=0.5):

    ## read in fasta file 2025-09-10-decoys-contam-UP000008827.fas
    fasta_obj = fasta.read("data/PRIDE/PXD005187/2025-09-11-decoys-contam-UP000005640.fas")
    fasta_lentghs = {}
    for description, sequence in fasta_obj:
        prefix, accession, rest = description.split('|', 2)
        if 'rev_' in prefix:
            continue
        fasta_lentghs[accession] = len(sequence)
    # print(fasta_lentghs)


    keep_df = []
    for name, name_df in df.groupby('ProteinName'):
        uniprot = name
        # print(name, uniprot)
        # print(name_df)

        seq = []
        start_end = name_df[['Protein.Start', 'Protein.End']].values
        # print(start_end)
        for s, e in start_end:
            seq.append(np.arange(s, e+1))
        seq = np.hstack(seq).flatten()
        seq = np.unique(seq)
        seq.sort()
        # print(f'Unique peptide positions for {uniprot}:', seq)

        if uniprot not in fasta_lentghs:
            # print(f'No match found in FASTA file for {uniprot}')
            quit()
        length = fasta_lentghs[uniprot]
        # print('Length of canonical sequence:', length)
        coverage = len(seq) / length * 100
        # print(f'Coverage for {uniprot}: {coverage:.2f}%')
        
        if coverage >= threshold * 100:
            name_df['Coverage'] = coverage
            keep_df.append(name_df)
        # else:
        #     print(f'Excluding {uniprot} due to low coverage: {coverage:.2f}%')

    if keep_df:
        keep_df = pd.concat(keep_df)
    else:
        keep_df = pd.DataFrame()  # Return empty DataFrame if no proteins pass the filter
    return keep_df
############################################################################

############################################################################
# 1) Load MSstats-ready table from IonQuant
df = pd.read_csv("data/Human_msstats.csv")
df = df[df["IsotopeLabelType"].isin(["H","L"])]
print(df)

# 1b) Check coverage
df = coverage(df, threshold=0.25)
print(df)

############################################################################

############################################################################
# 2) Sum to Protein × Run × Label, then pivot to H/L
summ = (df.groupby(["ProteinName","Run","IsotopeLabelType"], as_index=False)["Intensity"]
          .sum())
wide = (summ.pivot(index=["ProteinName","Run"],
                   columns="IsotopeLabelType",
                   values="Intensity")
            .reset_index())
print(wide)
############################################################################

############################################################################
# 3) Guard against zeros/Infs and compute log2(H/L)
wide = wide.replace([np.inf, -np.inf], np.nan)
zero_counts = []
for name, name_df in wide.groupby('ProteinName'):
    num_H_zeros = (name_df[['H']] == 0).sum().sum()
    num_L_zeros = (name_df[['L']] == 0).sum().sum()
    total_zeros = num_H_zeros + num_L_zeros
    zero_counts.append({'ProteinName': name, 'num_H_zeros': num_H_zeros, 'num_L_zeros': num_L_zeros, 'total_zeros': total_zeros})
zero_counts_df = pd.DataFrame(zero_counts)

LOD = 500000

H_zeros_loc = np.where(wide["H"] == 0)[0]
gaussian = np.random.normal(loc=LOD, scale=0.3, size=len(H_zeros_loc))
wide.loc[H_zeros_loc, "H"] = gaussian

L_zeros_loc = np.where(wide["L"] == 0)[0]
gaussian = np.random.normal(loc=LOD, scale=0.3, size=len(L_zeros_loc))
wide.loc[L_zeros_loc, "L"] = gaussian

wide = wide.dropna(subset=["H","L"])
wide = wide[(wide["H"] > 0) & (wide["L"] > 0)]
wide["log2HL"] = np.log2(wide["H"] / wide["L"])
print(wide)
############################################################################

############################################################################
# 4) save this intermediate protein level abundance file
wide.to_csv("data/Human_protein_abundance.csv", index=False)
print("Wrote data/Human_protein_abundance.csv")
############################################################################

############################################################################
# 5) Per-protein one-sample t-test of log2(H/L) vs 0
def ttest_stats(x):
    x = pd.Series(x).dropna()
    n = x.size
    mean = x.mean() if n else np.nan
    p = np.nan
    if n >= 2:
        _, p = ttest_1samp(x, 0.0)
    return pd.Series({"n_runs": n, "mean_log2HL": mean, "p_value": p})

# NOTE: the unstack() is the key change vs your output
res = (wide.groupby("ProteinName")["log2HL"]
          .apply(ttest_stats)
          .unstack())
print(res)
############################################################################

############################################################################
# 6) Benjamini–Hochberg FDR
from scipy.stats import false_discovery_control
p = res["p_value"]
res["q_value"] = false_discovery_control(p.values)
print(res)
############################################################################

############################################################################
# 7) Sort by q_value and print
out = res.reset_index().sort_values("q_value", na_position="last")
out = out.reset_index(drop=True)
sig = out[out['q_value'] < 0.05]
sig = sig.reset_index(drop=True)
print(f"Number of significant proteins (q < 0.05): {sig.shape[0]} out of {out.shape[0]}")
print(sig.to_string())
############################################################################

############################################################################
# 8) Save
out = res.reset_index().sort_values("p_value", na_position="last")
out.to_csv("data/Human_protein_abundance_log2HPLP_stats.csv", index=False)
print("Wrote data/Human_protein_abundance_log2HPLP_stats.csv")
############################################################################

############################################################################
# 9) plot qvalue vs mean_log2HP/LP
import matplotlib.pyplot as plt
# plt.scatter(out['mean_log2HP/LP'], -np.log10(out['q_value']))
plt.scatter(out['mean_log2HL'], -np.log10(out['p_value']))
plt.xlabel('mean_log2HL')
# plt.xlim(-6, 6)
plt.ylabel('-log10(q_value)')
plt.title('Volcano plot of Protein log2(H/L) vs -log10(q_value)')
plt.axhline(-np.log10(0.05), color='red', linestyle='--')
plt.savefig('data/Human_volcano_plot.png')
print("Wrote data/Human_volcano_plot.png")
############################################################################