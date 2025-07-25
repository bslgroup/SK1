#from Hope
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.stats import t
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from scipy.stats import ks_2samp
from scipy import stats
from matplotlib.patches import Patch
#Below is to have all of the data from the simulations
rep1_data = pd.read_csv("replicate1.txt", delim_whitespace=True)
rep2_data = pd.read_csv("replicate2.txt", delim_whitespace=True)
rep3_data = pd.read_csv("replicate3.txt", delim_whitespace=True)
rep4_data = pd.read_csv("replicate4.txt", delim_whitespace=True)
rep5_data = pd.read_csv("replicate5.txt", delim_whitespace=True)
rep6_data = pd.read_csv("replicate6.txt", delim_whitespace=True)
rep7_data = pd.read_csv("replicate7.txt", delim_whitespace=True)
rep8_data = pd.read_csv("replicate8.txt", delim_whitespace=True)
rep9_data = pd.read_csv("replicate9.txt", delim_whitespace=True)
rep10_data = pd.read_csv("replicate10.txt", delim_whitespace=True)
#Below is to exclue the first 100ns of the simulations
rep1_data_cut = rep1_data[rep1_data['frames'] > 50]
rep2_data_cut = rep2_data[rep2_data['frames'] > 50]
rep3_data_cut = rep3_data[rep3_data['frames'] > 50]
rep4_data_cut = rep4_data[rep4_data['frames'] > 100]
rep5_data_cut = rep5_data[rep5_data['frames'] > 100]
rep6_data_cut = rep6_data[rep6_data['frames'] > 100]
rep7_data_cut = rep7_data[rep7_data['frames'] > 100]
rep8_data_cut = rep8_data[rep8_data['frames'] > 100]
rep9_data_cut = rep9_data[rep9_data['frames'] > 100]
rep10_data_cut = rep10_data[rep10_data['frames'] > 100]

angstrom_symbol= "\u212B"

sns.set_style("whitegrid")
sns.set_context("poster")
#fig, axes = plt.subplots(2, 2, figsize=(12,12), share=True)
fig, ax = plt.subplots()

plt.xticks(fontsize=16)
colors = ['green'] * 3 + ['royalblue'] * 7
#sns.violinplot(y="rmsd2", data = del_data, ax=ax)
sns.violinplot(data=[rep1_data[["rmsd2"]],rep2_data[["rmsd2"]],rep3_data[["rmsd2"]],rep4_data[["rmsd2"]],rep5_data[["rmsd2"]],rep6_data[["rmsd2"]],rep7_data[["rmsd2"]],rep8_data[["rmsd2"]],rep9_data[["rmsd2"]],rep10_data[["rmsd2"]]], ax=ax, palette=colors, cut=0)
plt.xticks(ticks=np.arange(10), labels=np.arange(1,11))
plt.xlabel('Replicate Number', fontsize=16)
plt.ylabel('RMSD (' + angstrom_symbol + ')', fontsize=16)
plt.title('U1', fontsize=20)

legend_elements = [Patch(facecolor='green', edgecolor='black', label='4 µs'),Patch(facecolor='royalblue', edgecolor='black', label='500 ns')]
plt.legend(handles=legend_elements, title="Run Time", title_fontsize=12, loc='upper right', fontsize=12, bbox_to_anchor=(1.015, 1.07))

plt.savefig("7n6ureplicates.png", dpi=300, bbox_inches='tight')
plt.show()

#Plot 2

fig, ax = plt.subplots()

plt.xticks(fontsize=16)

colors = ['green'] * 3 + ['royalblue'] * 7
sns.violinplot(data=[rep1_data_cut[["rmsd2"]],rep2_data_cut[["rmsd2"]],rep3_data_cut[["rmsd2"]],rep4_data_cut[["rmsd2"]],rep5_data_cut[["rmsd2"]],rep6_data_cut[["rmsd2"]],rep7_data_cut[["rmsd2"]],rep8_data_cut[["rmsd2"]],rep9_data_cut[["rmsd2"]],rep10_data_cut[["rmsd2"]]], ax=ax, palette=colors, cut=0)
plt.xticks(ticks=np.arange(10), labels=np.arange(1,11))
plt.xlabel('Replicate Number', fontsize=16)
plt.ylabel('RMSD (' + angstrom_symbol + ')', fontsize=16)
plt.title('U1 Excluding the First 100 ns', fontsize=20)

legend_elements = [Patch(facecolor='green', edgecolor='black', label='4 µs'),Patch(facecolor='royalblue', edgecolor='black', label='500 ns')]
plt.legend(handles=legend_elements, title="Run Time", title_fontsize=12, loc='upper right', fontsize=12, bbox_to_anchor=(1.015, 1.03))

plt.savefig("7n6ureplicates_cut.png", dpi=300, bbox_inches='tight')
plt.show()


#sns.violinplot(x="Solubility", y="denovo_ddg_avg", order=["Soluble", "Insoluble"], data = del_data, ax=axes[0,0])
#sns.swarmplot(x="Solubility", y="denovo_ddg_avg", data = del_data, ax=axes[0,0], color="black", order=["Soluble", "Insoluble"])
#axes[0,0].set_title("DeNovo", size=20)
#axes[0,0].set_ylabel("$\Delta$$\Delta$G (REU)\n", size=20)
#axes[0,0].set_xlabel("", size=20)
#axes[0,0].tick_params(labelsize=18)
#axes[0,0].set_ylim([-20, 50])

#add_stat_annotation(axes[0,0], data=del_data, x="Solubility", y="denovo_ddg_avg", order=["Soluble", "Insoluble"], box_pairs=[("Soluble", "Insoluble")], test='Mann-Whitney', text_format='star', loc='inside',line_offset_to_box=0.2)


#sns.violinplot(x="Solubility", y="seghyb_ddg_avg", data = del_data, ax=axes[0,1], order=["Soluble", "Insoluble"])
#sns.swarmplot(x="Solubility", y="seghyb_ddg_avg", data = del_data, ax=axes[0,1], color="black", order=["Soluble", "Insoluble"])
#axes[0,1].set_title("Hybridize", size=20)
#axes[0,1].set_ylabel("", size=20)
#axes[0,1].set_xlabel("", size=20)
#axes[0,1].tick_params(labelsize=18)
#axes[0,1].set_ylim([-20, 50])

#add_stat_annotation(axes[0,1], data=del_data, x="Solubility", y="seghyb_ddg_avg", order=["Soluble", "Insoluble"], box_pairs=[("Soluble", "Insoluble")], test='Mann-Whitney', text_format='star', loc='inside',line_offset_to_box=0.15)


#sns.violinplot(x="Solubility", y="relax_ddg_avg", data = del_data, ax=axes[1,0], order=["Soluble", "Insoluble"])
#sns.swarmplot(x="Solubility", y="relax_ddg_avg", data = del_data, ax=axes[1,0], color="black", order=["Soluble", "Insoluble"])
#axes[1,0].set_title("Relax", size=20)
#axes[1,0].set_ylabel("$\Delta$$\Delta$G (REU)\n", size=20)
#axes[1,0].set_xlabel("", size=20)
#axes[1,0].set_xticklabels(["Soluble", "Insoluble"], size=20)
#axes[1,0].tick_params(labelsize=18)
#axes[1,0].set_ylim([-20, 50])

#add_stat_annotation(axes[1,0], data=del_data, x="Solubility", y="relax_ddg_avg", order=["Soluble", "Insoluble"], box_pairs=[("Soluble", "Insoluble")], test='Mann-Whitney', text_format='star', loc='inside',line_offset_to_box=0.15)


#sns.violinplot(x="Solubility", y="af_ddg_avg", data = del_data, ax=axes[1,1], order=["Soluble", "Insoluble"])
#sns.swarmplot(x="Solubility", y="af_ddg_avg", data = del_data, ax=axes[1,1], color="black", order=["Soluble", "Insoluble"])
#axes[1,1].set_title("AlphaFold", size=20)
#axes[1,1].set_ylabel("", size=20)
#axes[1,1].set_xlabel('', size=20)
#axes[1,1].set_xticklabels(["Soluble", "Insoluble"], size=20)
#axes[1,1].tick_params(labelsize=18)
#axes[1,1].set_ylim([-20, 50])

#add_stat_annotation(axes[1,1], data=del_data, x="Solubility", y="af_ddg_avg", order=["Soluble", "Insoluble"], box_pairs=[("Soluble", "Insoluble")], test='Mann-Whitney', text_format='star', loc='inside',line_offset_to_box=0.15)


#plt.savefig("all_2violin_af_avg_mannwhit_bigfont.png", dpi=300, bbox_inches='tight')
