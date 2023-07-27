import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import glob

# Plots generated 5/30/23


def plot_sample(path, ax, min_coverage=20, size=3, alpha=.2, spines=False):
    x = pd.read_table(path, index_col=0, header=None)
    x['coverage'] = x[4] + x[5]
    x = x[x['coverage'] > min_coverage]
    plt.scatter(x[1], x[3], s=size, alpha=alpha, c=x[3], cmap='jet',vmin=20, vmax=80,lw=0)
    plt.ylim([-10,110])
    ax.spines['top'].set_visible(spines)
    ax.spines['right'].set_visible(spines)   
    ax.set_xticks([])
    ax.set_yticks([0,50,100])
    plt.xlim([0, np.max(x[1])])
    return x


ctl = '/Volumes/de_novo/Bisulfite/Akata/coverage/C3.C3_1_bismark_bt2_pe.deduplicated.bismark.cov.EBV.cov'
igg = '/Volumes/de_novo/Bisulfite/Akata/coverage/Ig1.Ig1_1_bismark_bt2_pe.deduplicated.bismark.cov.EBV.cov'

# Whole genome
fig = plt.figure(figsize=(5,2))

ax = plt.subplot(211)
plot_sample(ctl, ax, 10)

ax2 = plt.subplot(212)
plot_sample(igg, ax2, 10)

plt.savefig('/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/whole_genome_latent_lytic_methylation.jpeg',dpi=1000)


def plot_sample_partial(path, ax, min_coverage=20, size=3, alpha=.2, spines=False, left=0, right=500000):
    x = pd.read_table(path, index_col=0, header=None)
    x['coverage'] = x[4] + x[5]
    x = x[x['coverage'] > min_coverage]
    x = x[(x[1]>=left) & (x[1] <= right)]
    plt.scatter(x[1], x[3], s=size, alpha=alpha, c=x[3], cmap='jet',vmin=20, vmax=80,lw=0)
    plt.ylim([-10,110])
    ax.spines['top'].set_visible(spines)
    ax.spines['right'].set_visible(spines)   
    ax.set_xticks([])
    ax.set_yticks([0,50,100])
    plt.xlim([left, right])




# Qp-EBNA1 Zoom
fig = plt.figure(figsize=(2,2))
ax = plt.subplot(111)
plot_sample_partial(ctl, ax, 10,40,1,True, 113104, 113404)

plt.savefig('/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/Qp-EBNA1_bisulfite_zoomin.svg',dpi=1000)
plt.savefig('/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/Qp-EBNA1_bisulfite_zoomin.jpeg',dpi=1000)

# Cp-EBNA1 Zoom
fig = plt.figure(figsize=(2,2))
ax = plt.subplot(111)
plot_sample_partial(ctl, ax, 10,40,1,True,74300,74600)
ax.set_xlim([74300, 74600])
plt.savefig('/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/Cp-EBNA1_bisulfite_zoomin.svg',dpi=1000)
plt.savefig('/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/Cp-EBNA1_bisulfite_zoomin.jpeg',dpi=1000)

## Get methylated CpG fraction by gene

bed = pd.read_table('/Users/nate/Documents/Genomes/Beds/chrEBV_Akata_inverted_refined_genes_plus_features_annotation_cleaned.bed', index_col=0, skiprows=1, header=None)
bed = bed.set_index(3)

words = ['Repeat', 'splice', 'Bind', 'miR', 'repeat', 'DRlef', 'DRr']
throw = []
for i in bed.index:
    for word in words:
        if word in i:
            throw.append(i)
bed = bed.loc[[i for i in bed.index if i not in throw]]
pos = bed[bed[5] == '+']
neg = bed[bed[5] == '-']

latent_files = glob.glob('/Volumes/de_novo/Bisulfite/Akata/coverage/C*EBV.cov')
promoter_length = 200
min_counts = 10
df = pd.DataFrame(index=bed.index)
num_events = pd.DataFrame(index=bed.index)
min_sites = 1

for path in latent_files:
    cov = pd.read_table(path, index_col=0, header=None)
    cov['coverage'] = cov[4] + cov[5]
    cov = cov[cov['coverage'] >= min_counts]
    d = {}
    number_of_events = {}
    for i in pos.index:
        start = pos.loc[i, 1] - promoter_length
        stop = pos.loc[i, 1]
        subset = cov[(cov[1] >= start) & (cov[1] <= stop)][3]
        d[i] = np.mean(subset)
        number_of_events[i] = len(subset.index)
    
    for i in neg.index:
        start = neg.loc[i, 2] 
        stop = neg.loc[i, 2] + promoter_length      
        subset = cov[(cov[2] >= start) & (cov[2] <= stop)][3]
        d[i] = np.mean(subset)
        number_of_events[i] = len(subset.index)


    df[path.split('/')[-1].split('.')[0]] = [d[i] for i in df.index]
    num_events[path.split('/')[-1].split('.')[0]] = [number_of_events[i] for i in num_events.index]

num_events = num_events[num_events>min_sites].dropna()
df = df.loc[num_events.index]
df['mean'] = np.mean(df,1)
df = df.sort_values('mean')

fig = plt.figure(figsize=(12,4))
ax = plt.subplot()
for x in range(len(df.index)):
    plt.scatter([x, x+.01, x+.02, x+.03],df.iloc[x,:4],s=2,alpha=.7, c='k', ec='k')
    plt.bar(x, df.iloc[x,4], color='k', alpha=.5, lw=.3, ec='k')


plt.ylim([0,100])
plt.yticks([0,50,100])
plt.xlim([-.5, len(df.index)])
plt.savefig("/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure4/bisulfite_seq_plots/ebv_gene_methyl_percentage.svg")