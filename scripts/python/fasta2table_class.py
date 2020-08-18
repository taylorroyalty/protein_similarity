# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 07:00:47 2020

@author: Taylor Royalty
"""

from Bio import SeqIO
from re import split as r_split

class fasta_anno_analysis():
    
    # init method or constructor
    def __init__(self,fasta,db_source):
        self.fasta = SeqIO.parse(open(fasta),"fasta")
        self.db_source= db_source
        
    def swiss_fasta2table(self):

        def swiss_anno_seq(fasta_inter):
            annotation, sequence= fasta_inter.description, str(fasta_inter.seq)
            annotation=(annotation.split(" OS=")[0]).split(" ", 1)[1]
            return annotation,sequence

        self.anno_seq=[swiss_anno_seq(entry) for entry in self.fasta]

    def ncbi_fasta2table(self):

        def ncbi_anno_seq(fasta_inter):
            annotation, sequence= fasta_inter.description, str(fasta_inter.seq)
            annotation=(r_split(" \[.*\]$",annotation)[0]).split(" ", 1)[1]
            return annotation,sequence

        return [ncbi_anno_seq(entry) for entry in self.fasta]





