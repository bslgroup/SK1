import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === PARAMETERS ===
rep_dirs = ["../1/", "../2/", "../3/"]
input_file = "lipid_contact_CA_summary.txt"
output_png = "lipid_contact_heatmap_colorful_avg.png"

# === HARDCODED number of frames per replicate ===
rep_frame_counts = {
    "../1/": 573,
    "../2/": 590,
    "../3/": 620
}

# === LOAD & NORMALIZE EACH REPLICATE ===
dfs = []

for rep in rep_dirs:
    path = os.path.join(rep, input_file)
    total_frames = rep_frame_counts[rep]
    
    df = pd.read_csv(path, sep='\t')
    lipid_types = df.columns[1:]
    
    for lipid in lipid_types:
        df[lipid] = df[lipid] / total_frames * 100  # normalize to % per frame count
    
    dfs.append(df)

# === COMPUTE AVERAGE ACROSS REPLICATES ===
df_avg = dfs[0].copy()
for lipid in lipid_types:
    df_avg[lipid] = sum(df[lipid] for df in dfs) / len(dfs)

# === PLOT HEATMAP ===
df_heatmap = df_avg.set_index("Residue")[lipid_types]

plt.figure(figsize=(10, 4))
sns.heatmap(
    df_heatmap,
    cmap="crest",  # or try 'plasma', 'magma', 'rocket'
    linewidths=0,
    cbar_kws={'label': 'Avg % of Simulation Time'},
    xticklabels=True,
    yticklabels=True
)

plt.title("Residue–Lipid Contact (Avg % of Simulation Time, 3 Reps)")
plt.xlabel("Lipid")
plt.ylabel("Residue")
plt.tight_layout()
plt.savefig(output_png, dpi=300)
print(f"✅ Heatmap saved to {output_png}")
plt.show()
