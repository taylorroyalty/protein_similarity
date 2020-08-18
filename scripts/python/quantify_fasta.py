# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 07:00:47 2020

@author: Taylor Royalty
"""

#%%

import csv

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
        annotation, sequence= seq_rec.description, str(seq_rec.seq)
        annotation=(annotation.split(" OS=")[0]).split(" ", 1)[1]
        return [annotation,sequence]
    
    fasta_iter = SeqIO.parse(open(fasta),"fasta") #convert fasta file to fastaiterator object
    anno_seq_list=[swiss_anno_seq(entry) for entry in fasta_iter] #return a list of sequences and annotations
    
    if write == True:
        with open(filename,'w',newline='') as f:
            wr=csv.writer(f,delimiter=delim)
            wr.writerow(['annotation','sequence'])
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

def count_sequence_motifs(aa_motifs,sequences):
    