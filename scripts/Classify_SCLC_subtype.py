import os
import pandas as pd
import re
import numpy as np
import glob
from copy import deepcopy
import sys
import phenograph
import scanpy as sc
import itertools

out_dir = '/data/peer/chanj3/HTA.SCLC.010920/out.publication.010920/'
os.makedirs(out_dir, exist_ok=True)

### Can be downloaded at https://data.humantumoratlas.org/
adata = sc.read_h5ad(out_dir + 'adata.SCLC.010920.h5ad')

'''
INITIAL SUBTYPE SCORING
'''

SCLC_subtype_dir='/data/peer/chanj3/ref/bulk.SCLC/' ### CHANGE TO DIRECTORY HOLDING BULK DEG FOR SCLC SUBTYPE

ct_dict = {}
for i in ['SCLC-A','SCLC-N','SCLC-P','SCLC-Y']:
    SCLC_subtype_file = SCLC_subtype_dir + 'George_etal.limma.DEG.%s.090219.csv' % i

    deg_rest = pd.read_csv(SCLC_subtype_file,sep='\t')
    ind = deg_rest.logFC > 0

    ct_dict[i] = deg_rest.loc[ind, 'adj.P.Val'].sort_values().index[:100]


adata_tmp = deepcopy(adata)

sc.pp.scale(adata_tmp)

for i in ct_dict.keys():
    sc.tl.score_genes(adata_tmp, gene_list = set(ct_dict[i]),score_name = 'limma_' + i, use_raw = False) #Note that scores can change marginally  based on randomly selected reference genes

'''
GET TRAINING LABELS FOR CLASSIFICATION
'''

SCLCtype_score = adata_tmp.obs.loc[:,['limma_SCLC-A','limma_SCLC-N','limma_SCLC-P']]


extrema = pd.DataFrame([[10,0,0],[0,10,0],[0,0,10]],
                          index = ['SCLC-A', 'SCLC-N', 'SCLC-P'],
                          columns = ['limma_SCLC-A','limma_SCLC-N','limma_SCLC-P'])

# Correlation with single cells
from sklearn.metrics import pairwise_distances

# Distances
dists = pd.DataFrame(pairwise_distances(extrema,
                                        SCLCtype_score, metric='euclidean'),
                     index=extrema.index, columns=SCLCtype_score.index)

cell_types = dists.index

train_labels = {}
for cell_type in cell_types:
    train_labels[cell_type] = dists.loc[cell_type,:].sort_values().iloc[:100].index #Top 100 from score_genes

print((pd.Index(itertools.chain(*[train_labels[i] for i in train_labels.keys()])).value_counts() > 1).any()) #Ensure that there are no overlaps for top training labels across subtypes

'''
MARKOV ABSORPTION CLASSIFICATION USING PHENOGRAPH.CLASSIFY
'''

num_dcs = 10 #Can adjust how many DCs to retain for Phenograph classification
dm_ev = pd.DataFrame(adata.obsm['X_diffmap'], index = adata.obs.index).loc[:,:num_dcs]
train = np.empty((len(cell_types),),dtype=object)

for c,cell_type in enumerate(cell_types):
    labels = train_labels[cell_type]
    train[c] = dm_ev.loc[labels,:]
test = dm_ev

SCLCtype_pval = phenograph.classify(train, test, k=30, metric='euclidean')


adata.obs = adata.obs.drop(adata.obs.columns[adata.obs.columns.str.contains('pval_')], axis=1)

adata.obs = pd.concat([adata.obs,
                       pd.DataFrame(SCLCtype_pval[1], index = dm_ev.index,
                                    columns = [re.sub('limma_','pval_',i) for i in SCLCtype_score.columns])],
                      axis=1)

adata.obs.loc[:,'SCLC_subtype'] = adata.obs.loc[:,['pval_SCLC-A', 'pval_SCLC-N', 'pval_SCLC-P']].idxmax(axis=1).str.replace('pval_','')
