# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:49:35 2020

@author: Taylor Royalty
"""

from sklearn.cluster import OPTICS
from joblib import Parallel, delayed
# from sklearn.decomposition import PCA
from sklearn.manifold import MDS
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import math
import os
import warnings


dist_dir='data/lev_distances/'
write_path='data/mds_cluster_dataframes/'
n=1
pnt=27

def parallel_OPTICS_cluster_sensitivity(filename,dist_dir,pnt,write_path):
    #suppress runtime warning for OPTICS algorithm. dividing by small number 
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    file=np.loadtxt(open(dist_dir+filename,'rt').readlines()[:-1],delimiter=",")
    file=np.loadtxt(dist_dir+filename,delimiter=',')[:-1,:] 
    id_list=file[-1,:]
    file=file[:-1,:]
    

    #perform multidimensional scaling on the levenshtein distance matrix for 2D visualization
    mds_obj=MDS(n_components=2, dissimilarity='precomputed', random_state=1)
    mds_comp=mds_obj.fit_transform(file)

    # i=0
    #create stoarge array for results
    # df_tmp=pd.DataFrame(data=np.zeros(shape=(10,3)),columns=["annotation","clusters","minpoints"])
    # df_tmp["minpoints"]=minpoint_range
    # df_tmp["annotation"]=filename
    print(filename)
    # for pnt in minpoint_range:    
    optic_obj=OPTICS(min_samples=pnt,algorithm='ball_tree')
    optic_clust=optic_obj.fit(file)
    clust_df=pd.DataFrame(data=mds_comp,columns=["Component_1","Component_2"])
    clust_df['Cluster']=optic_clust.labels_
    clust_df['swiss_id']=id_list
        # df_tmp.iloc[i,1]=len(np.unique(optic_clust.labels_[optic_clust.labels_>-1]))
        # i+=1
    
    clust_df.to_csv(write_path+filename, header=True, index=False)


Parallel(n_jobs=n)(delayed(parallel_OPTICS_cluster_sensitivity)(filename,dist_dir,pnt,write_path) for filename in os.listdir(dist_dir))

    
    # n_clust_list.append([filename,n_clust_tmp,minpoint_range])
    # with open('data/cluster_analysis','w',newline='') as f:
    #     wr=csv.writer(f,delimiter=',')
    #     # wr.writerow(['annotation','c'])
    #     wr.writerows(n_clust_list)
        
        
    # plt.scatter(minpoint_range,n_clust)

    # clust_df.to_csv(write_path+filename,sep=",",index=False)

