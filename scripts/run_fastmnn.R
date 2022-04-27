suppressMessages(library(Matrix))
suppressMessages(library(batchelor))
suppressMessages(library(BiocParallel))
ncores = 16
SnowParam(workers = ncores)
suppressMessages(library(irlba))
library(dplyr)

knn <- 30 #as.integer(args[1])
num_hvg <- 5000 # as.integer(args[2])

dir = '/data/peer/chanj3/HTA.combined.010920/out.combined.010920/' #CHANGE FOR YOUR WORKING DIRECTORY

#NUMBER OF HVGs PRE-CALCULATED
hvg_fn = sprintf('%shvg.%d.txt', num_hvg)
hvg = readLines(hvg_fn)


file = paste0(dir, 'norm_dbtf.combined.barcodes.010920.csv')
bc = read.csv(file,header=F,col.names = 'Cell', stringsAsFactors = F)
batch = gsub('_[0-9]+$','',bc$V1)

file = paste0(dir, 'norm_dbtf.combined.genes.010920.csv')
g = read.csv(file,header=F,col.names = 'Gene', stringsAsFactors = F)

file = paste0(dir, 'norm_dbtf.combined.010920.mtx') #READ IN SPARSE MATRIX OF NORMALIZED COUNTS. Median normalized used in SCLC atlas
counts0 = as.matrix(readMM(file))

rownames(counts0) = bc$Cell
colnames(counts0) = g$Gene
counts0 = t(counts0)

#READ IN METADATA (AVAILABLE ON HTAN DATA PORTAL)
obs_df = read.table('/home/chanj3/data/HTA.combined.010920/out.mnnc.010920/obs.combined.010920.txt', sep ='\t', header=T)
obs_df = obs_df %>% dplyr::select(batch, patient, histo, tissue)
obs_df2 = obs_df[!duplicated(obs_df$batch),] %>% dplyr::arrange(batch)

num_tissue = sapply(sort(unique(obs_df2$patient)), function(x) length(unique(obs_df2[obs_df2$patient==x, 'tissue'])))
num_cell = table(obs_df$batch)

obs_df$num_tissue = num_tissue[obs_df$patient]
obs_df$num_cell = num_cell[obs_df$batch]
                    
meta = obs_df %>% dplyr::rename(sample=batch)                
meta = meta[bc$Cell,]

order_df = meta[!duplicated(meta$sample), ]
order_df$histo = gsub('normal','LUAD',order_df$histo) #RELABEL NORMAL TO LUAD SINCE THEY ARE NORMAL LUNG ADJACENT TO LUAD

#HIERARCHICAL MERGING STRATEGY: MERGE FIRST BY PATIENT THEN BY HISTOLOGY. MERGE PATIENTS WITH MULTIPLE TISSUE SITES FIRST. WITHIN PATIENT, MERGE SAMPLES WITH HIGHEST NUMBER OF CELLS FIRST
order_df = order_df[order(order_df$histo, order_df$num_tissue, order_df$patient, order_df$ncells, decreasing = TRUE),]

batch_order = sort(unique(order_df$sample))

merge.order = list()
for (i in unique(order_df$patient)) {
    merge.order[[i]] = which(order_df$patient==i)
}

merge.order2 = list()
for (i in unique(order_df$histo)) {
    tmp = order_df$patient[order_df$histo==i]
    merge.order2[[i]] = merge.order[tmp[!duplicated(tmp)]]
}

counts = counts0[hvg,]

batch = gsub('_[0-9]+$','',colnames(counts))

counts_list = lapply(split(seq_along(batch), batch), function(m, ind) m[,ind], m=counts)[order(unique(batch))]

counts_list = counts_list[batch_order]

out = fastMNN(counts_list, batch = batch_order, k=knn, merge.order = merge.order, cos.norm=F, d = 50, get.variance=T)

correct = reducedDims(out)$corrected
correct = correct[match(colnames(counts), rownames(correct)),]
rotation = rowData(out)$rotation

correct.per_gene <- tcrossprod(correct,rotation)

odir= paste0(dir, 'fastMNN_output/')

corrpc_ofile = paste0('fastMNN.KNN_',as.character(knn),'.HVG_',as.character(num_hvg),'.corrected_pc.csv')
correxp_ofile = paste0('fastMNN.KNN_',as.character(knn),'.HVG_',as.character(num_hvg),'.corrected_exp.csv')

corrpc_ofile = paste0(odir,corrpc_ofile)
correxp_ofile = paste0(odir,correxp_ofile)

write.table(correct, file=corrpc_ofile,sep="\t")
write.table(correct.per_gene, file=correxp_ofile,sep="\t")
