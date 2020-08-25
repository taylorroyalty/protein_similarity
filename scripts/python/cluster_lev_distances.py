# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:49:35 2020

@author: Taylor Royalty
"""

from sklearn.cluster import OPTICS
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from matplotlib.pyplot import plot as plt
import numpy as np

file=np.loadtxt('data/lev_distances/ATP_synthase_subunit_delta.csv',delimiter=',', dtype=float) 

optic_obj=OPTICS(min_samples=20)
pca_obj=PCA(n_components=2)
mds_obj=MDS(n_components=2, dissimilarity='precomputed', random_state=1)

optic_clust=optic_obj.fit(file)
pca_comp=pca_obj.fit(file)
mds_comp=mds_obj.fit(file)

plt.scatter(mds_comp[:, 0], mds_comp[:, 1], **colorize)