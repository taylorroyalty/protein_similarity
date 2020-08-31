# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:49:35 2020

@author: Taylor Royalty
"""

from sklearn.cluster import OPTICS
from joblib import Parallel, delayed
# from sklearn.decomposition import PCA
# from sklearn.manifold import MDS
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import os
import warnings


dist_dir='data/lev_distances/'
write_path='data/cluster_analysis/'
n=1

def parallel_OPTICS_cluster_sensitivity(filename,dist_dir,write_path):
    #suppress runtime warning for OPTICS algorithm. dividing by small number 
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    file=np.loadtxt(open(dist_dir+filename,'rt').readlines()[:-1],delimiter=",")
    file=file[:-1,:]
    
    #determine number of samples
    n_samp=file.shape[1]
    #create a vector of evenly lognormally-spaced minpoints to test for minimum size of the clster
    minpoint_range=2**np.linspace(math.log2(5),math.log2(n_samp/10),20)
    minpoint_range=minpoint_range.astype(int)

    i=0
    #create stoarge array for results
    df_tmp=pd.DataFrame(data=np.zeros(shape=(10,3)),columns=["annotation","clusters","minpoints"])
    df_tmp["minpoints"]=minpoint_range
    df_tmp["annotation"]=filename
    print(filename)
    for pnt in minpoint_range:    
        optic_obj=OPTICS(min_samples=pnt,algorithm='ball_tree')
        optic_clust=optic_obj.fit(file)
        df_tmp.iloc[i,1]=len(np.unique(optic_clust.labels_[optic_clust.labels_>-1]))
        i+=1
    
    df_tmp.to_csv(write_path+filename, header=True, index=False)

Parallel(n_jobs=n)(delayed(parallel_OPTICS_cluster_sensitivity)(filename,dist_dir,write_path) for filename in os.listdir(dist_dir))

    
    # n_clust_list.append([filename,n_clust_tmp,minpoint_range])
    # with open('data/cluster_analysis','w',newline='') as f:
    #     wr=csv.writer(f,delimiter=',')
    #     # wr.writerow(['annotation','c'])
    #     wr.writerows(n_clust_list)
        
        
    # plt.scatter(minpoint_range,n_clust)
    # clust_df=pd.DataFrame(data=mds_comp,columns=["Component_1","Component_2"])
    # clust_df['Cluster']=optic_clust.labels_
    # clust_df.to_csv(write_path+filename,sep=",",index=False)

