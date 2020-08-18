# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 07:00:47 2020

@author: Taylor Royalty
"""

from Bio import SeqIO
from re import split as r_split
import csv 
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
