# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:49:35 2020

@author: Taylor Royalty
"""

from sklearn.cluster import OPTICS
from joblib import Parallel, delayed
# from sklearn.manifold import MDS
from sklearn.manifold import TSNE

import numpy as np
import pandas as pd
import os
import warnings


dist_dir='data/lev_distances/'
write_path='data/mds_cluster_dataframes/'
seq_path='data/swiss_n100.tsv'
n=1
pnt=27

swiss100=pd.read_csv(seq_path,sep='\t')

def parallel_OPTICS_cluster_sensitivity(filename,dist_dir,pnt,write_path):
    #suppress runtime warning for OPTICS algorithm. dividing by small number 
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    file=np.loadtxt(open(dist_dir+filename,'rt').readlines()[:-1],delimiter=",")
    id_list=open(dist_dir+filename,'rt').readlines()[-1].split("\n")[0].split(",")
    

    #perform multidimensional scaling on the levenshtein distance matrix for 2D visualization
    tsne = TSNE(n_components=2, random_state=0)
    tsne_comp=tsne.fit_transform(file)
    # mds_obj=MDS(n_components=2, dissimilarity='precomputed', random_state=1)
    # mds_comp=mds_obj.fit_transform(file)

    print(filename)
    # perform OPTICS clustering 
    optic_obj=OPTICS(min_samples=pnt,algorithm='ball_tree')
    optic_clust=optic_obj.fit(file)
    
    #generate data.frame
    clust_df=pd.DataFrame(data=tsne_comp,columns=['Component_1','Component_2'])
    clust_df['Cluster']=optic_clust.labels_
    clust_df['id']=id_list
    
    #merge sequence data from 
    clust_df=pd.merge(clust_df,swiss100,on='id')

    
    clust_df.to_csv(write_path+filename, header=True, index=False)


Parallel(n_jobs=n)(delayed(parallel_OPTICS_cluster_sensitivity)(filename,dist_dir,pnt,write_path) for filename in os.listdir(dist_dir))

    
    # n_clust_list.append([filename,n_clust_tmp,minpoint_range])
    # with open('data/cluster_analysis','w',newline='') as f:
    #     wr=csv.writer(f,delimiter=',')
    #     # wr.writerow(['annotation','c'])
    #     wr.writerows(n_clust_list)
        
        
    # plt.scatter(minpoint_range,n_clust)

    # clust_df.to_csv(write_path+filename,sep=",",index=False)

