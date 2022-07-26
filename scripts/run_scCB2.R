suppressMessages(require("scCB2"))
suppressMessages(require("Matrix"))
args <- commandArgs(trailingOnly = TRUE)

# command line arguments
mtx_file <- args[1]
dir <- args[2]

ofile = gsub('_sparse_molecule_counts.mtx', '_realdrops.csv',mtx_file)

bc_file = gsub('molecule_counts.mtx','counts_barcodes.csv',mtx_file)
barcodes = read.csv(paste0(dir,bc_file),stringsAsFactors = F,
                        header=F,sep=',',row.names=1,as.is=T)
barcodes = as.character(barcodes$V2)

g_file = gsub('molecule_counts.mtx','counts_genes.csv',mtx_file)
genes = read.csv(paste0(dir,g_file),stringsAsFactors = F,
                      header=F,sep=',',row.names=1,as.is=T)
genes = as.character(genes$V2)

my.counts = readMM(paste0(dir,mtx_file))

rownames(my.counts) <- barcodes
colnames(my.counts) <- genes
my.counts = t(my.counts)

ind = !duplicated(genes)
my.counts = my.counts[ind,]

my.counts = as(my.counts,'dgCMatrix')

CBOut <- CB2FindCell(my.counts, FDR_threshold = 0.01,
    lower = 100, Ncores = 2, verbose = TRUE)
summary(CBOut)

RealCell <- GetCellMat(CBOut)
str(RealCell)

real_cells = colnames(RealCell)

write.table( data.frame(real_cells), paste0(dir,ofile), quote = F, row.names = F, col.names = F )
