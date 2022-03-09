import pandas as pd
import numpy as np

def chrName2Int(chr):
    if chr=="X":                #    encode chr X as 23
        return 23
    elif chr=="Y":              #    encode chr Y as 24
        return 24
    else:
        return int(chr)

def gene_cmp(a, b, geneD):
    if geneD[a][0] > geneD[b][0]:
        return 1
    elif geneD[a][0] < geneD[b][0]:
        return -1
    else:
        if geneD[a][1] > geneD[b][1]:
            return 1
        elif geneD[a][1] < geneD[b][1]:
            return -1
        else:
            return 0

def cmp_to_key(mycmp, geneD):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj, geneD) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj, geneD) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj, geneD) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj, geneD) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj, geneD) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj, geneD) != 0
    return K

def infer_cnv(genePosFile, chrFile, expdf, narounds=50, method='genes'):

    print("Number of around bins/genes to be taken in account: "+str(narounds))

def infer_cnv(genePosFile, chrFile, expdf, narounds=50, method='genes'):

    print("Number of around bins/genes to be taken in account: "+str(narounds))
    print("Moving average by "+method)

    geneD=dict()

    # Read gene position file
    with open(genePosFile, "r") as f:
        for line in f:
            arr=line.split("\t")
            if arr[1]!="MT":
                if (arr[0] not in geneD):
                    geneD[arr[0]]=[chrName2Int(arr[1]),int(arr[2]),int(arr[3])]
                else:
                    start=geneD[arr[0]][1]
                    end=geneD[arr[0]][2]
                    if abs(end-start)<abs(int(arr[2])-int(arr[3])):
                        geneD[arr[0]]=[chrName2Int(arr[1]),int(arr[2]),int(arr[3])]

    #    Read chromosome file
    chrL=dict()
    with open(chrFile, "r") as f:
        for line in f:
            arr=line.split("\t")
            if (arr[0].isnumeric()) or (arr[0]=="X") or (arr[0]=="Y"):
                chrL[chrName2Int(arr[0])]=int(arr[1])
                #print(arr[0]+" "+arr[1])
    #expdf=pd.read_csv(expFile,index_col=0)

    geneArr=expdf.columns

    print("Missing genes if any:")
    ct = 0
    for g in geneArr:
        if g not in geneD:
            #print(g)
            expdf = expdf.drop(g, 1)
            ct+=1
    print(ct)

    geneArr=list(expdf.columns)
    geneArr=sorted(list(geneArr),key=cmp_to_key(gene_cmp, geneD))
    expdf=expdf[geneArr]


    print("starting moving average to infer cnv")

    binGArr=[0]*24
    for i in range(len(geneArr)):
        binGArr[geneD[geneArr[i]][0]-1]+=1

    chrArr=[]
    for i in range(24):
        for j in range(binGArr[i]):
            chrArr.append(i+1)

    expdf.columns = [chrArr, expdf.columns]

    for i in range(1, 24):
        if i == 1:
            cnvdf = expdf[i].rolling(window=(2*narounds)+1, center=True, min_periods=0, axis=1).mean()
        else:
            cnvdf = pd.concat([cnvdf, expdf[i].rolling(window=(2*narounds)+1, center=True, min_periods=0, axis=1).mean()], axis=1)

    print("done moving average to infer cnv")
    print(len(chrArr))
    print(len(cnvdf.columns))
    cnvdf.columns=[chrArr[:len(cnvdf.columns)], cnvdf.columns]
    #cnvdf.to_csv(outputFile)

    allv=cnvdf.values.flatten()
    allv=allv[np.nonzero(allv)]

    print("* median non-zero values: "+str(np.percentile(allv,50)))
    print("* 5-percentile of non-zero values: "+str(np.percentile(allv,5)))
    print("* 95-percentile of non-zero values: "+str(np.percentile(allv,95)))

    return cnvdf
