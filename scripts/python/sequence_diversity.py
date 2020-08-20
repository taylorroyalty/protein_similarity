# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:24:15 2020

@author: Taylor
"""

#import packages and modules
import sys
import quantify_fasta as qf
import pandas as pd

#relative path to python scripts
sys.path.insert(1,'scripts/python')

#read in fasta file and convert to dataframe
swiss = qf.swiss_fasta2table('data/swiss_prot_08032020.fasta')
swiss = pd.DataFrame(swiss,columns=(['annotation','sequence']))

#count annotations in swiss prot
swiss_max = swiss.annotation.mode()
swiss_max_sequence=swiss[swiss.annotation == swiss_max[0]]

#define amino acid motifs based on sequence length
aa_motif=qf.unique_aa_motifs(3)


motif_count_df=qf.count_sequence_motifs(aa_motif,swiss_max_sequence['sequence'].tolist())

motif_sums=motif_count_df.groupby(by='motif').sum()
