import sys
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

def plotcnv(cnvArr, chrArr, title, norm=None, norm_vals=None, norm_inds=None,
            cmap=plt.cm.bwr, minv=None, maxv=None, save_name=None, collapse_to_zero=None,
            collapse_max_to=None, collapse_min_to=None, dendrogram=True, linkage_method='ward'):
    fig = plt.figure(figsize=(16,8))
    if dendrogram:
        ax1 = fig.add_axes([0.01,0.1,0.1,0.6])
    if not isinstance(cnvArr, pd.DataFrame):
        cnvdf=pd.DataFrame(cnvArr)
    else:
        cnvdf = cnvArr
    newcnvdf = cnvdf.copy()
    if norm == 'mean':
        if isinstance(norm_inds, np.ndarray):
            newcnvdf = newcnvdf.subtract(newcnvdf.loc[norm_inds].mean(axis=0),axis=1).as_matrix()
        elif isinstance(norm_vals, np.ndarray):
            newcnvdf = newcnvdf
        else:
            newcnvdf=newcnvdf.subtract(newcnvdf.mean(axis=0),axis=1).as_matrix()
    elif norm == 'mode':
        newcnvdf = newcnvdf.subtract(newcnvdf.mode(axis=0),axis=1).as_matrix()
    else:
        newcnvdf = newcnvdf.as_matrix()

    if dendrogram:
        Y = sch.linkage(newcnvdf, method=linkage_method)
        Z1 = sch.dendrogram(Y, truncate_mode='level', orientation='right')
        ax1.set_xticks([])
        ax1.set_yticks([])


    # identify maximum value as 95 percentile
    #allv=cnvdf.values.flatten()

    allv=np.ravel(newcnvdf)
    allv=allv[np.nonzero(allv)]
    if maxv == None and collapse_max_to == None:
        maxv=np.percentile(allv,99)
    if minv == None and collapse_min_to == None:
        minv=np.percentile(allv,1)
    if maxv == None:
        maxv = collapse_max_to
    if minv == None:
        minv = collapse_min_to

    if norm != 'none':
        if -minv>maxv:
            maxv=-minv
        if maxv>-minv:
            minv=-maxv

    if collapse_min_to != None:
        inds = np.where(newcnvdf < collapse_min_to)
        newcnvdf[inds] = np.full(len(inds[0]), collapse_min_to)
    if collapse_max_to != None:
        inds = np.where(newcnvdf > collapse_max_to)
        newcnvdf[inds] = np.full(len(inds[0]), collapse_max_to)
    if collapse_to_zero != None:
        inds = np.where((newcnvdf >= -collapse_to_zero) & (newcnvdf <= collapse_to_zero))
        newcnvdf[inds] = np.zeros(len(inds[0]))
        inds = np.where(newcnvdf > collapse_to_zero)
        newcnvdf[inds] = np.subtract(newcnvdf[inds], collapse_to_zero)
        inds = np.where(newcnvdf < -collapse_to_zero)
        newcnvdf[inds] = np.add(newcnvdf[inds], collapse_to_zero)

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.12,0.1,0.7,0.6])
    if dendrogram:
        idx1 = np.array(sch.leaves_list(Y))
        newcnvdf = newcnvdf[idx1,:]
    im = axmatrix.matshow(newcnvdf, aspect='auto', origin='lower',
                                                  cmap=cmap,vmin=minv,vmax=maxv)

    xticks=[]
    xcords=[]
    xticklabels=[]
    for i in range(0,len(chrArr)):
        if (i==len(chrArr)-1) or (chrArr[i]!=chrArr[i+1]):
            xcords.append(i)
            plt.axvline(x=i,color='black')
            if int(chrArr[i])==23:
                xticklabels.append('X')
            elif int(chrArr[i])==24:
                xticklabels.append('Y')
            else:
                xticklabels.append(chrArr[i])
    for i in range(len(xcords)):
        if i==0:
            xticks.append(xcords[i]/2.0)
        else:
            xticks.append(xcords[i-1]+(xcords[i]-xcords[i-1])/2.0)
    plt.xticks(xticks,xticklabels)
    axmatrix.set_yticks([])
    axmatrix.set_title(title)

    # Plot colorbar.
    axcolor = fig.add_axes([0.83,0.1,0.02,0.6])
    cbar = plt.colorbar(im, cax=axcolor)
    if collapse_to_zero != None:
        ticks = cbar.ax.yaxis.get_ticklabels()
        new_ticks = []
        for i in range(len(ticks)):
            curTick = ticks[i].get_text()
            if i < np.floor(len(ticks)/2):
                new_ticks.append(-1*(float(curTick[1:]) + collapse_to_zero))
            elif i > np.floor(len(ticks)/2):
                new_ticks.append(float(curTick) + collapse_to_zero)
            else:
                new_ticks.append(0.0)
        cbar.ax.set_yticklabels(new_ticks)
    # plt.show()

    if save_name != None:
        plt.savefig(save_name, bbox_inches='tight')
# plot 2 cnv arrays together with reference one refCnvArr on top

"""
def plot2cnv(cnvArr, refCnvArr, chrArr, title, refTitle, minv=0.0,maxv=0.5):
    cnvdf=pd.DataFrame(cnvArr)
    refcnvdf=pd.DataFrame(refCnvArr)
    #refCnvArr=refcnvdf.subtract(cnvdf.mean(axis=0),axis=1).as_matrix()
    #cnvArr=cnvdf.subtract(cnvdf.mean(axis=0),axis=1).as_matrix()
    refCnvArr=refcnvdf.subtract(refcnvdf.mean(axis=0),axis=1).as_matrix()
    cnvArr=cnvdf.subtract(cnvdf.mean(axis=0),axis=1).as_matrix()

    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_axes([0.01,0.1,0.1,0.6])
    Y = sch.linkage(cnvArr, method='ward')
    Z1 = sch.dendrogram(Y, orientation='right')
    Y2 = sch.linkage(refCnvArr, method='ward')
    Z2 = sch.dendrogram(Y2, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])



    allv=np.concatenate((np.ravel(cnvArr),np.ravel(refCnvArr)),axis=0)
    allv=allv[np.nonzero(allv)]
    maxv=np.percentile(allv,98)
    minv=np.percentile(allv,2)

    # computing x ticks
    xticks=[]
    xcords=[]
    xticklabels=[]
    for i in range(0,len(chrArr)):
        if (i==len(chrArr)-1) or (chrArr[i]!=chrArr[i+1]):
            xcords.append(i)
            #plt.axvline(x=i,color='black')
            if int(chrArr[i])==23:
                xticklabels.append('X')
            elif int(chrArr[i])==24:
                xticklabels.append('Y')
            else:
                xticklabels.append(chrArr[i])
    for i in range(len(xcords)):
        if i==0:
            xticks.append(xcords[i]/2.0)
        else:
            xticks.append(xcords[i-1]+(xcords[i]-xcords[i-1])/2.0)

    axmatrix = fig.add_axes([0.12,0.1,0.7,0.6])
    idx1 = Z1['leaves']
    cnvArr = cnvArr[idx1,:]
    im = axmatrix.matshow(cnvArr, aspect='auto', origin='lower', cmap=plt.cm.bwr,vmin=minv,vmax=maxv)
    for i in range(0,len(chrArr)):
        if (i==len(chrArr)-1) or (chrArr[i]!=chrArr[i+1]):
            plt.axvline(x=i,color='black')
    plt.xticks(xticks,xticklabels)
    axmatrix.set_yticks([])
    axmatrix.set_title(title)

    axrefmatrix = fig.add_axes([0.12,0.8,0.7,0.6])
    idx2 = Z2['leaves']
    refCnvArr = refCnvArr[idx2,:]
    im = axrefmatrix.matshow(refCnvArr, aspect='auto', origin='lower', cmap=plt.cm.bwr,vmin=minv,vmax=maxv)
    for i in range(0,len(chrArr)):
        if (i==len(chrArr)-1) or (chrArr[i]!=chrArr[i+1]):
            plt.axvline(x=i,color='black')
    plt.xticks(xticks,xticklabels)
    axrefmatrix.set_yticks([])
    axrefmatrix.set_title(refTitle)

    # Plot colorbar.
    axcolor = fig.add_axes([0.83,0.1,0.02,0.6])
    plt.colorbar(im, cax=axcolor)
    
"""
