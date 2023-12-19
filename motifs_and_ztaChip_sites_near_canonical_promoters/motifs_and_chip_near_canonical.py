import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
import pyBigWig


# awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$2+1,1}' Raji_Zta_ChIP_induced_pooledReps_summits.hg38.bed > raji_chip.bg
# sort -k1,1 -k2,2n raji_chip.bg > raji_chip.sorted.bg
# bedGraphToBigWig raji_chip.sorted.bg ~/Documents/Genomes/star/hg38/chrNameLength.txt raji_chip.sorted.bw



upstream_bases = 5500
downstream_bases = 5500
length = upstream_bases + downstream_bases

canonical_over3tpm_bed_path = '/Users/nate/dnovo/beds/TSS_from_mutu_ctl2_over3tpm.bed'
canonical_0tpm_bed_path = '/Users/nate/dnovo/beds/mutu_canonical_TSS_0tpm.bed'
denovo_bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.noblacklist.bed'
zta_chip_bw_path = '/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure3/starting_files/raji_chip.sorted.bw'
ap1_motif_bw_path = '/Users/nate/dnovo/beds/TGAGTCA.bw'
meZta_motif_bw_path = '/Users/nate/dnovo/beds/TGAGCGA.bw'
other_ap1_motif_bw_path = '/Users/nate/dnovo/beds/TGACTCA.bw'


def import_bed(path, col_sort=False):
    regions = []
    with open(path) as bed_handle:
        for line in bed_handle:
            regions.append(line.strip('\n').split('\t'))
    random.shuffle(regions)
    if col_sort:
        regions.sort(key=lambda x: float(x[col_sort]))
    return regions


def bin(array_2d, bin_size):
    l = []
    for i in range(0, len(array_2d), bin_size):
        l.append(np.sum(array_2d[i:i+bin_size]))
    return l
    
    
def moving_average(a, n=10):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def extract_coverage(bw_path, regions):
    
    bw = pyBigWig.open(bw_path)
    m = np.zeros([len(regions), length])
    for ind, region in enumerate(regions):
        start = int(float(region[1]))
        stop = int(float(region[2])) 
        strand = region[5]
        middle = (start + stop) // 2 
        start = middle - upstream_bases
        stop = middle + downstream_bases
        vals = bw.values(region[0], start, stop)
        vals = np.abs(vals)
        vals = np.nan_to_num(np.array(vals), nan=0)
        if strand == '-':
            vals = np.flip(vals)
        m[ind] +=  vals
    print("done")
    return m 


def main(bed_path, bw_path, binsize=250, m_avg=4):
    regions = import_bed(bed_path)
    cov = extract_coverage(bw_path, regions=regions)
    cov_sum = np.sum(cov, 0)
    binned = bin(cov_sum, binsize)
    return binned
    #return moving_average(binned, m_avg)

def plot_lines(lines=[], names=[], binsize=250):
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot()
    xrange = [ i + binsize/2 for i in range(-upstream_bases, downstream_bases, binsize)]
    for line, name in zip(lines, names):
        line = line / np.mean(line[0:5])
        plt.plot(xrange, line, label=name)
        plt.fill_between(xrange, line, 0.9*np.min(line), alpha=.3)

    ax.set_xlim([-5000, 5000])
    ax.set_ylim([0.9*np.min(line),1.1*np.max(line)])
    ax.set_xticks([-5000, 0, 5000])
    plt.axvline(0, c='k', ls='--',alpha=.5)
    return fig



for bedpath, bedname in zip([canonical_over3tpm_bed_path, canonical_0tpm_bed_path, denovo_bed_path], ['over3tpm', '0tpm', 'dnovo']):
    for bwpath, bwname in zip([zta_chip_bw_path, ap1_motif_bw_path,meZta_motif_bw_path,other_ap1_motif_bw_path ], ['zta_chip', 'TGAGTCA', 'TGAGCGA', 'TGACTCA']):
        line = main(bedpath, bwpath)
        fig = plot_lines([line],['name'])
        plt.savefig(f'{bedname}_{bwname}.svg')
        plt.close()
