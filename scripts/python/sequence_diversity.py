# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:24:15 2020

@author: Peng
"""
#compile function generating sequences

import sys
sys.path.insert(1,'scripts/python')

from fasta2table import swiss_fasta2table

swiss=swiss_fasta2table('data/swiss_prot_08032020.fasta')