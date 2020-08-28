# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 13:31:43 2020

@author: Peng
"""

from Bio import SeqIO
import csv

def tara_fasta2table(fasta,write=False,filename='output.tsv',delim='\t'):
# =============================================================================
#--fasta is a text file containing sequence data in fasta format. This handles
#   multiline fasta. This has only been varied to work on ncbi fasta files
#--write is a boolean which controls whether the output list is written to a file
#   default is False
#--filename is the name of the file to be written. The default is output.tsv
# =============================================================================

    def anno_seq(seq_rec):
# =============================================================================
#--seq_rec is a sequence record from a fastaiterator object.
#   the sequence record is decomposed into the sequence and annotation. 
#   the annotation is parsed from the annotation.
# =============================================================================
        unique_id, sequence=seq_rec.id,  str(seq_rec.seq)
        return [unique_id, sequence]
  
        
    fastaiter = SeqIO.parse(open(fasta),"fasta")#convert fasta file to fastaiterator object
    anno_seq_list = [anno_seq(entry) for entry in fastaiter]#return a list of sequences and annotations    
    
    if write == True:
        with open(filename,'w',newline='') as f:
            wr=csv.writer(f,delimiter=delim)
            wr.writerow(['id','sequence'])
            wr.writerows(anno_seq_list)
            
    return anno_seq_list
