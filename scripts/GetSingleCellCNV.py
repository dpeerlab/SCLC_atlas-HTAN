# load modules
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict, Counter, OrderedDict
from imp import reload
from itertools import combinations
from copy import deepcopy
import matplotlib.pyplot as plt

from plotCnvArr import plotcnv
from inferCNVFast import infer_cnv

import seaborn as sns
import re
import scanpy as sc

mode='NSCLC_epithelial'
out_dir = sys.argv[1]
adata = sys.argv[2]
cnv_file= 'cnv_df.NSCLC_epithelial.010920.h5'
plot_file = 'cnv.hc.NSCLC_epithelial.010920.png'

sns.set_style('white')

norm_unfiltered_df = pd.DataFrame(adata.layers['norm_df'], index=adata.obs_names, columns = adata.var_names)

data_df = np.log2(norm_unfiltered_df+0.1)

batch = [re.sub('_[0-9]+$','',i) for i in data_df.index]
batch = pd.Series(batch, index = data_df.index)

genePosFile = '/home/chanj3/data/HTA.NSCLC_epithelial.plasticity.070219/notebooks/infercnv/genePos.wo_MT.txt'
chrFile = '/home/chanj3/data/HTA.NSCLC_epithelial.plasticity.070219/notebooks/infercnv/chrNameLength.txt'

cnv_df = infer_cnv(genePosFile,chrFile,data_df, narounds=100)

##### REGEX '_N_|_LN_minus' OBTAINS NORMAL/WILDTYPE
diploid_mean = cnv_df.loc[cnv_df.index.str.contains('_N_|_LN_minus'),:].mean(axis=0)
diploid_std = cnv_df.loc[cnv_df.index.str.contains('_N_|_LN_minus'),:].std(axis=0)

cnv_df2 = cnv_df.subtract(diploid_mean, axis=1).div(diploid_std, axis=1)


allv=np.ravel(cnv_df2)
allv=allv[np.nonzero(allv)]
maxv=np.percentile(allv,99)
minv=np.percentile(allv,1)
print(maxv, minv)

cnv_df2 = cnv_df2.sub(cnv_df2.median(axis=1), axis=0)
wt_std = cnv_df2.loc[cnv_df2.index.str.contains('_N_|_LN_minus'),:].std(axis=0)

from copy import deepcopy
tmp = deepcopy(cnv_df2.values)
tmp[np.where(np.abs(cnv_df2).lt(1.5*wt_std, axis=1))] = 0
cnv_df3 = pd.DataFrame(tmp, index=cnv_df2.index, columns=cnv_df2.columns)

store = pd.HDFStore(out_dir + cnv_file)
store['cnv_df'] = cnv_df3.astype(np.float32)
store.close()

batch = [re.sub('_[0-9]+$','',i) for i in cnv_df3.index]
batch = pd.Series(batch, index = cnv_df3.index)

pal_tn = sns.color_palette('Set1', 2)
lut_tn = dict(zip([False, True], pal_tn))

pal_tp = sns.color_palette('Set1', len(set(batch)))
lut_tp = dict(zip(sorted(list(set(batch)), key = lambda x: '_N_|_LN_minus' in x), pal_tp))

pal_chr = sns.color_palette('deep', 23)
lut_chr = dict(zip(list(range(1,24)), pal_chr))

row_colors1 = [lut_tp[ re.sub('_[0-9]+$','',i) ] for i in cnv_df3.index]
row_colors2 = [lut_tn[ '_N_' in i or 'LN_minus' in i ] for i in cnv_df3.index]
col_colors = [lut_chr[i] for i in cnv_df3.columns.get_level_values(0)]

def normalize(data):
    means = np.ravel(data.mean(axis=1))
    std = np.ravel(data.std(axis=1))
    ret = ((data.T - means) / std).T
    return ret

g = sns.clustermap(cnv_df3, col_cluster = False, center = 0, vmin=-4, vmax=4, method='ward',
                   cmap = plt.cm.seismic,
                   col_colors = col_colors,
                   row_colors = [row_colors1,row_colors2])
ax = g.ax_heatmap
ax.set_xlabel("")
ax.set_xticklabels("")

fname = out_dir + plot_file
g.savefig(fname, dpi=300)

