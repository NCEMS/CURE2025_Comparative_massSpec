import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
from pyteomics import fasta

############################################################################
def coverage(df, threshold=0.5):

    ## read in fasta file 2025-09-10-decoys-contam-UP000008827.fas
    fasta_obj = fasta.read("data/PRIDE/PXD023343/2025-09-10-decoys-contam-UP000008827.fas")

    keep_df = []
    for name, name_df in df.groupby('ProteinName'):
        uniprot = name.split('|')[1]
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

        length = None
        for description, sequence in fasta_obj:
            if uniprot in description:
                # print('Found match in FASTA file')
                # print('Description:', description)
                # print('Sequence:', sequence)
                length = len(sequence)
                # print('Length of sequence:', length)
                break
        if length is None:
            # print(f'No match found in FASTA file for {uniprot}')
            continue
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
df = pd.read_csv("data/Soybean_msstats.csv")
df = df[df["IsotopeLabelType"].isin(["H","L"])]

# 1a) Map condition to HP or LP
condition = ["HP" if "HP" in c else "LP" for c in df["Condition"]]
df["Condition"] = condition
run = [c[-1] for c in df["Run"]]
df["Run"] = run
print(df)

# 1b) Check coverage
df = coverage(df, threshold=0.25)
print(df)
############################################################################

############################################################################
# 2) Sum to Protein × Run × Label, then pivot to H/L
summ = (df.groupby(["ProteinName","Run","Condition"], as_index=False)["Intensity"].sum())
wide = (summ.pivot(index=["ProteinName","Run"],
                   columns="Condition",
                   values="Intensity").reset_index())
############################################################################

############################################################################
# 3) Guard against zeros/Infs and compute log2(H/L)
wide = wide.replace([np.inf, -np.inf], np.nan)
zero_counts = []
for name, name_df in wide.groupby('ProteinName'):
    num_HP_zeros = (name_df[['HP']] == 0).sum().sum()
    num_LP_zeros = (name_df[['LP']] == 0).sum().sum()
    total_zeros = num_HP_zeros + num_LP_zeros
    zero_counts.append({'ProteinName': name, 'num_HP_zeros': num_HP_zeros, 'num_LP_zeros': num_LP_zeros, 'total_zeros': total_zeros})
zero_counts_df = pd.DataFrame(zero_counts)

LOD = 500000

HP_zeros_loc = np.where(wide["HP"] == 0)[0]
gaussian = np.random.normal(loc=LOD, scale=0.3, size=len(HP_zeros_loc))
wide.loc[HP_zeros_loc, "HP"] = gaussian

LP_zeros_loc = np.where(wide["LP"] == 0)[0]
gaussian = np.random.normal(loc=LOD, scale=0.3, size=len(LP_zeros_loc))
wide.loc[LP_zeros_loc, "LP"] = gaussian

wide = wide.dropna(subset=["HP","LP"])
wide = wide[(wide["HP"] > 0) & (wide["LP"] > 0)]
wide["log2HP/LP"] = np.log2(wide["HP"] / wide["LP"])
print(wide)

uniprots = [n.split('|')[1] for n in wide['ProteinName']]
wide['ProteinName'] = uniprots
print(wide)
############################################################################

############################################################################
# 4) save this intermediate protein level abundance file
wide.to_csv("data/Soybean_protein_abundance.csv", index=False)
print("Wrote data/Soybean_protein_abundance.csv")
############################################################################

############################################################################
# 5) t-test on log2(H/L) for each protein
def ttest_stats(x):
    x = pd.Series(x).dropna()
    n = x.size
    mean = x.mean() if n else np.nan
    # standard error
    stderr = x.std(ddof=1) / np.sqrt(n) if n >= 2 else np.nan
    p = np.nan
    if n >= 2:
        _, p = ttest_1samp(x, 0.0)
    return pd.Series({"n_runs": n, "mean_log2HP/LP": mean, "stderr_log2HP/LP": stderr, "p_value": p})


# NOTE: the unstack() is the key change vs your output
res = (wide.groupby("ProteinName")["log2HP/LP"]
          .apply(ttest_stats)
          .unstack())
print(res.head())
# merge zero counts
res = res.merge(zero_counts_df, on='ProteinName', how='left')
print(res)
print(res['total_zeros'].value_counts())
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
out.to_csv("data/Soybean_protein_abundance_log2HPLP_stats.csv", index=False)
print("Wrote data/Soybean_protein_abundance_log2HPLP_stats.csv")
############################################################################

############################################################################
# 9) plot qvalue vs mean_log2HP/LP
import matplotlib.pyplot as plt
# plt.scatter(out['mean_log2HP/LP'], -np.log10(out['q_value']))
plt.scatter(out['mean_log2HP/LP'], -np.log10(out['p_value']))
plt.xlabel('mean_log2HP/LP')
# plt.xlim(-6, 6)
plt.ylabel('-log10(q_value)')
plt.title('Volcano plot of Protein log2(HP/LP) vs -log10(q_value)')
plt.axhline(-np.log10(0.05), color='red', linestyle='--')
plt.savefig('data/Soybean_volcano_plot.png')
print("Wrote data/Soybean_volcano_plot.png")
############################################################################