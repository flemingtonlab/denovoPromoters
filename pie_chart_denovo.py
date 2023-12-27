import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Primary analysis should be bedtools intersect
# Then make bed file with only the gene start sites and gene start site+1 for bedtools closest analysis

### WITHIN 1.5 KB CLOSEST ANALYSIS ###

# bedtools closest -s -D a -t first -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b /Users/nate/hg38_gene_coords.bed > denovo_closest_canonical_sense.bed
# bedtools closest -S -D a -t first -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b /Users/nate/hg38_gene_coords.bed  > denovo_closest_canonical_antisense.bed

# bedtools closest -s -D a -t first -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b /Users/nate/hg38_gene_TSSs.bed > denovo_closest_canonical_sense.bed
# bedtools closest -S -D a -t first -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b /Users/nate/hg38_gene_TSSs.bed > denovo_closest_canonical_antisense.bed

sense = pd.read_table("/Users/nate/denovo_closest_canonical_sense.bed", index_col=0, header=None)
antisense = pd.read_table("/Users/nate/denovo_closest_canonical_antisense.bed", index_col=0, header=None)
sense = sense.set_index(3)
antisense = antisense.set_index(3)
new = pd.DataFrame(index=sense.index)
new['sense'] = sense[12]
new['antisense'] = antisense[12]
new['upstream_sense_1500'] = new['sense'].map(lambda x:1 if 25 < x < 1500 else 0)
new['antisense_1500'] = new['antisense'].map(lambda x:1 if 25 < x < 1500 else 0)
new['both1500'] = np.sum(new[['upstream_sense_1500','antisense_1500']],1)

### GENE BODY INTERSECT ###

# awk '$3=="gene" {print $1,$4,$5,$14,"100",$7}' /Users/nate/Documents/Genomes/GTF/gencode.v44.basic.annotation.gtf |sed 's/"\(.*\).*";/\1/'|sed 's/ /  /g'>hg38_gene_coords.bed
# bedtools intersect -wa -b /Users/nate/hg38_gene_coords.bed  -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed  -u -s>denovo_intersect_with_genebody.sense.bed
# bedtools intersect -wa -b /Users/nate/hg38_gene_coords.bed  -a /Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed  -u -S>denovo_intersect_with_genebody.antisense.bed

dn = pd.read_table("/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed", index_col=3, header=None)
sense = pd.read_table("denovo_intersect_with_genebody.sense.bed", header=None, index_col=0)
antisense = pd.read_table("denovo_intersect_with_genebody.antisense.bed", header=None, index_col=0)
sense = sense.set_index(3)
antisense = antisense.set_index(3)
new2 = pd.DataFrame(index=dn.index)
new2['genebody_sense'] = sense[4]
new2['genebody_antisense'] = antisense[4]
new2[new2 > 1] = 1
new2 = new2.fillna(0)
new2['genebody_both'] = np.sum(new2[['genebody_sense','genebody_antisense']],1)

new[new2.columns] = new2




# Intergenic -> outside genebody and further than 1500bases away from gene
# Sense, upstream -> within 1500 bases of TSS
# Antisense, gene body 
# Sense, gene body

intergenic = new[(new['both1500']==0) & (new['genebody_both']==0)].shape[0]  # Intergenic: 6789
gb_sense = new[(new['genebody_antisense']==0) & (new['genebody_sense']==1)].shape[0]  # Gene body sense: 11200
gb_antisense = new[(new['genebody_antisense']==1) & (new['genebody_sense']==0)].shape[0]  # Gene body antisense: 8375
gb_both = new[(new['genebody_antisense']==1) & (new['genebody_sense']==1)].shape[0]  # Gene body both sense and antisense: 2060
upstream1500 = new[(new['genebody_both']==0) & (new['upstream_sense_1500']==1)].shape[0]  # Not in gene body, upstream 1500 and sense to gene: 1192 

labels = ["Sense", "Antisense", "Sense and antisense", "Sense (From -1500 to -25)", "Intergenic (>1500 bases away)"]
data = [gb_sense, gb_antisense, gb_both, upstream1500, intergenic]

fig = plt.figure()
ax = plt.subplot()
plt.bar(range(len(labels)), data, .7)
plt.savefig('bargraph_denovo_promoter_location.svg')