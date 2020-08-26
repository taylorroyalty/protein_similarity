# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:49:35 2020

@author: Taylor Royalty
"""

from sklearn.cluster import OPTICS
# from sklearn.decomposition import PCA
from sklearn.manifold import MDS
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import os
import csv

dist_dir='data/lev_distances/'
write_path='data/mds_cluster_dataframes/'

n_clust_list=[]
for filename in os.listdir(dist_dir):
    
    file=np.loadtxt(dist_dir+filename,delimiter=',', dtype=float) 
    
    #determine number of samples
    n_samp=file.shape[1]
    #create a vector of evenly lognormally-spaced minpoints to test for minimum size of the clster
    minpoint_range=10**np.linspace(math.log10(5),math.log10(n_samp/10),10)
    minpoint_range=minpoint_range.astype(int)
    
    #perform multidimensional scaling on the levenshtein distance matrix for 2D visualization
    mds_obj=MDS(n_components=2, dissimilarity='precomputed', random_state=1)
    mds_comp=mds_obj.fit_transform(file)
    
    i=0
    #create stoarge array for results
    n_clust_tmp=np.zeros(shape=(10,),dtype=int)
    for pnt in minpoint_range:    
        optic_obj=OPTICS(min_samples=pnt,algorithm='ball_tree',n_jobs=1)
        optic_clust=optic_obj.fit(file)
        n_clust_tmp[0,i]=len(np.unique(optic_clust.labels_[optic_clust.labels_>-1]))
        i+=1
    
    n_clust_list.append([filename,n_clust_tmp,minpoint_range])
    with open('data/cluster_analysis','w',newline='') as f:
        wr=csv.writer(f,delimiter=',')
        # wr.writerow(['annotation','c'])
        wr.writerows(n_clust_list)
    # plt.scatter(minpoint_range,n_clust)
    # clust_df=pd.DataFrame(data=mds_comp,columns=["Component_1","Component_2"])
    # clust_df['Cluster']=optic_clust.labels_
    # clust_df.to_csv(write_path+filename,sep=",",index=False)

