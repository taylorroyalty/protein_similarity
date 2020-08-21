# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 07:00:47 2020

@author: Taylor Royalty
"""

#%%

import csv
import pandas as pd
import numpy as np
import Levenshtein as lev #C libraries; much faster than python script
#https://codereview.stackexchange.com/questions/217065/calculate-levenshtein-distance-between-two-strings-in-python
#https://pypi.org/project/python-Levenshtein/
from itertools import combinations_with_replacement
from Bio import SeqIO
from re import split as r_split
 
#%%    
def swiss_fasta2table(fasta,write=False,filename='output.tsv',delim='\t'):    
# =============================================================================
#--fasta is a text file containing sequence data in fasta format. This handles
#   multiline fasta. This has only been varied to work on swiss prot fasta files
#--write is a boolean which controls whether the output list is written to a file
#   default is False
#--filename is the name of the file to be written. The default is output.tsv
# =============================================================================

    def swiss_anno_seq(seq_rec):
# =============================================================================
#--seq_rec is a sequence record from a fastaiterator object.
#   the sequence record is decomposed into the sequence and annotation. 
#   the annotation is parsed from the annotation.
# =============================================================================
        unique_id, annotation, sequence=seq_rec.id, seq_rec.description, str(seq_rec.seq)
        annotation=(annotation.split(" OS=")[0]).split(" ", 1)[1]
        return [unique_id, annotation,sequence]
    
    fasta_iter = SeqIO.parse(open(fasta),"fasta") #convert fasta file to fastaiterator object
    anno_seq_list=[swiss_anno_seq(entry) for entry in fasta_iter] #return a list of sequences and annotations
    
    if write == True:
        with open(filename,'w',newline='') as f:
            wr=csv.writer(f,delimiter=delim)
            wr.writerow(['id','annotation','sequence'])
            wr.writerows(anno_seq_list)
            
    return anno_seq_list

#%%

def ncbi_fasta2table(fasta,write=False,filename='output.tsv',delim='\t'):
# =============================================================================
#--fasta is a text file containing sequence data in fasta format. This handles
#   multiline fasta. This has only been varied to work on ncbi fasta files
#--write is a boolean which controls whether the output list is written to a file
#   default is False
#--filename is the name of the file to be written. The default is output.tsv
# =============================================================================

    def ncbi_anno_seq(seq_rec):
# =============================================================================
#--seq_rec is a sequence record from a fastaiterator object.
#   the sequence record is decomposed into the sequence and annotation. 
#   the annotation is parsed from the annotation.
# =============================================================================
        annotation, sequence= seq_rec.description, str(seq_rec.seq)
        annotation=(r_split(" \[.*\]$",annotation)[0]).split(" ", 1)[1]
        return annotation,sequence
        
    fastaiter = SeqIO.parse(open(fasta),"fasta")#convert fasta file to fastaiterator object
    anno_seq_list = [ncbi_anno_seq(entry) for entry in fastaiter]#return a list of sequences and annotations    
    
    if write == True:
        with open(filename,'w',newline='') as f:
            wr=csv.writer(f,delimiter=delim)
            wr.writerow(['annotation','sequence'])
            wr.writerows(anno_seq_list)
            
    return anno_seq_list

#%%

def unique_aa_motifs(n):
# =============================================================================
# generates a list containing all unique amino acid motifs.
# --n is an integer which defines length of motifs
# =============================================================================
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    
    aa_motifs = list(combinations_with_replacement(aa,n))
    aa_motifs_forward = [''.join(i) for i in aa_motifs] 
    aa_motifs_backward = [i[::-1] for i in aa_motifs_forward] #considers reverse order
    aa_motif_all = aa_motifs_forward+aa_motifs_backward #combine forward backwards for all combinations
    
    aa_motif_all =list(set(aa_motif_all))#set removes duplicate sequences, e.g., AA, YY, etc 
    
    return aa_motif_all

#%%

def count_sequence_motifs(aa_motif,sequences):
# =============================================================================
# a function for counting motifs in sequences. returns a Dataframe (columns=[motif, sequence, counts])
# --aa_motif is a list where elements are amino acid motifs being tabulated in seqeunces
# --sequences is a list containing amino acid sequences where motifs are being tabulated
# =============================================================================
    def CountOccurrences(main_string, substring):
# =============================================================================
# a subfunction for counts the occurence of substrings in string.
# --the string variable is the string to be counted for substring
# --substring is the string pattern being searched for in main_string
# =============================================================================
        count = 0
        start = main_string.find(substring, 0) 
        while start < len(main_string): 
            pos = main_string.find(substring, start)   
            if pos != -1: 
                start = pos + 1
                count += 1
            else: 
                break
        return count
    
    count_matrix=[[CountOccurrences(seq,motif) for seq in sequences]+[motif] for motif in aa_motif]
    count_df=pd.DataFrame(data=count_matrix,columns=(sequences+['motif']))
    return pd.melt(count_df,id_vars='motif',var_name='sequence',value_name='counts')

    
#%%
def iterative_levenshtein(s, t):
    #this function was written from https://www.python-course.eu/levenshtein_distance.php
    """ 
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings 
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein 
        distance between the first i characters of s and the 
        first j characters of t
    """

    rows = len(s)+1
    cols = len(t)+1
    dist = [[0 for x in range(cols)] for x in range(rows)]

    # source prefixes can be transformed into empty strings 
    # by deletions:
    for i in range(1, rows):
        dist[i][0] = i

    # target prefixes can be created from an empty source string
    # by inserting the characters
    for i in range(1, cols):
        dist[0][i] = i
        
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0
            else:
                cost = 1
            dist[row][col] = min(dist[row-1][col] + 1,      # deletion
                                 dist[row][col-1] + 1,      # insertion
                                 dist[row-1][col-1] + cost) # substitution

    # for r in range(rows):
    #     print(dist[r])
    
 
    return dist[row][col]

#%%
def parallel_lev_dist(seq_df,anno,write_path='data/lev_distances/'):
    tmp_df=seq_df[seq_df["annotation"] == anno].reset_index(drop=True)
    seq_n=len(tmp_df)
    lev_dist=np.zeros((seq_n,seq_n),dtype=int)
    # lev_dist=pd.DataFrame(data=np.zero(shape=(seq_n,seq_n)),columns=df["sequence"],index=df["sequence"])
    for i in range(seq_n):
        seq1=tmp_df["sequence"][i]
        for j in range(i,seq_n):
            seq2=tmp_df["sequence"][j]
            lev_dist[i,j]=lev.distance(seq1,seq2)
    
    #make similarity matrix symmetric
    lev_dist = lev_dist + lev_dist.T - np.diag(np.diag(lev_dist))
    
    #format results as dataframe 
    lev_dist_df=pd.DataFrame(data=lev_dist,columns=tmp_df["id"])
    lev_dist_df['id']=tmp_df["id"]
    lev_dist_df['annotation']=anno
    lev_dist_df=pd.melt(lev_dist_df,id_vars='id',var_name='id2',value_name='lev_dist')
    
    #write dataframe
    anno_correct=anno.replace('/','_')
    write_path=write_path+anno_correct+'.tsv'
    write_path=write_path.replace(' ','_')
    print(write_path)
    lev_dist_df.to_csv(write_path, header=True, index=False, sep='\t')
