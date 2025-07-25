import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and align replicates
res_order = ["LYS27", "LYS29", "ARG186", "LEU194", "PHE197", "LEU198"]
replicate_data = []

for rep in [1, 2, 3]:
    df = pd.read_csv(f"../{rep}/lipid_contacts_percentage.txt", sep="\t")
    df = df.set_index("Residue").loc[res_order].reset_index()
    replicate_data.append(df["Contact_Percentage"])

# Combine into single dataframe
df_combined = pd.DataFrame(replicate_data).T
df_combined.columns = [f"rep{r}" for r in [1, 2, 3]]
df_combined["Residue"] = res_order
df_combined["Mean"] = df_combined[["rep1", "rep2", "rep3"]].mean(axis=1)
df_combined["SEM"] = df_combined[["rep1", "rep2", "rep3"]].sem(axis=1)

# Plot settings
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "axes.spines.bottom": False,
})

fig, ax = plt.subplots(figsize=(9, 5))
ax.set_facecolor("white")

# Main color
bar_color = "#9B59B6"
sem_color = "#555555"  # darker gray for better visibility

# Draw bars with refined SEM styling
bars = ax.barh(df_combined["Residue"], df_combined["Mean"],
               xerr=df_combined["SEM"],
               color=bar_color, height=0.6, edgecolor='none',
               error_kw=dict(ecolor=sem_color, lw=1.2, capsize=4, alpha=0.8),
               zorder=3)

# Grid and axes
ax.invert_yaxis()
ax.xaxis.grid(True, linestyle='--', linewidth=0.6, alpha=0.5)
ax.set_xlim(0, 100)  # Set X-axis to 0–100%
ax.set_xlabel("Contact Percentage (%)", fontsize=14)
ax.set_title("Lipid Contact Frequency of Key Protein Residues (Mean ± SEM)", fontsize=16, weight='bold')

plt.tight_layout()
plt.savefig("lipid_contact_sem_final.pdf", dpi=300)
plt.show()
