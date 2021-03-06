# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 08:36:59 2020

@author: Taylor
"""

import sys
#relative path to python scripts
sys.path.insert(1,'scripts/python')

import quantify_fasta as qf
import pandas as pd

from joblib import Parallel, delayed


#number of cores to use
n = 1

#read in fasta file and convert to dataframe
swiss = qf.swiss_fasta2table('data/swiss_prot_08032020.fasta')
swiss = pd.DataFrame(swiss,columns=(['id','annotation','sequence']))

#generate all unique annotations
swiss_anno_500=swiss.groupby("annotation").filter(lambda x: len(x)>100).reset_index(drop=True)
swiss_anno_500.to_csv('data/swiss_n100.tsv',sep='\t',index=False)
swiss_anno_uniq=swiss_anno_500.annotation.unique()

#parallelize levenshstein distance calculation for fasta files
Parallel(n_jobs=n)(delayed(qf.parallel_lev_dist)(swiss_anno_500,anno,write_path='data/lev_distances_null/',null=True) for anno in swiss_anno_uniq)

# for anno in swiss_anno_uniq[0:2]:
#     swiss_tmp=swiss_anno_500[swiss_anno_500["annotation"] == anno]
#     seq_n=len(swiss_tmp)
#     lev_dist=np.zeros((seq_n,seq_n),dtype=int)
#     for i in range(seq_n):
#         seq1=swiss_tmp["sequence"][i]
#         for j in range(i+1,seq_n):
#             seq2=swiss_tmp["sequence"][j]
#             lev_dist[i,j]=qf.iterative_levenshtein(seq1,seq2)

