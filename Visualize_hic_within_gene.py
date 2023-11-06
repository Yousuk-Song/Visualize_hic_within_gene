#!/home/eaststar0/miniconda3/bin/python

import sys
import os
from tadlib.visualize.heatmaps import *

mcool = sys.argv[1]
loop = sys.argv[2]
tad = sys.argv[3]
input_gene = sys.argv[4]
res=10000

pos_file = open(f'{sys.path[0]}/gencode_protein_coding_gene.hg38.tsv', 'r')
n = 0
D = {}
for line in pos_file:
        if n == 0:
                n = 1
                continue
        cols = line.rstrip().split()

        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        strand = cols[3]
        gene = cols[4]
        D[gene] = [strand, chrom, start, end]
        if input_gene == gene:
                break
pos_file.close()

plot_window = 1000000
[strand, chrom, start, end] = D[input_gene]
print('gene: ', input_gene, chrom, start, end)
start = max(start - plot_window, 0)
end = end + plot_window

print(f'bin: {chrom}:{start}-{end}')

if norm == 'Y':
        vis = Triangle(f'{mcool}::/resolutions/{res}', chrom, start, end)
elif norm == 'N':
        vis = Triangle(f'{mcool}::/resolutions/{res}', chrom, start, end, None)
vis.matrix_plot()

vis.plot_TAD(tad, linewidth=1.5)
#vis.plot_loops(loop)

if not os.path.exists(gene):
        os.system(f'mkdir {gene}')

vis.outfig(gene + '/' + mcool.split('/')[-1].replace('.mcool', f'.vis.{gene}.png'))
