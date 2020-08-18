from itertools import combinations_with_replacement
import pandas as pd
import numpy as np

r=2
aa='ACDEFGHIKLMNPQRSTVWY'

sequences=pd.read_table('data/cytochrome_sequences',squeeze=True)

aa_motifs=list(combinations_with_replacement(aa,r))
aa_motifs_forward=np.array([''.join(i) for i in aa_motifs],dtype='str')
aa_motifs_backward=np.array([i[::-1] for i in aa_motifs_forward],dtype='str')
duplicate_aa=list(map(len,map(set,aa_motifs)))
index_duplicate=[idx for idx, val in enumerate(duplicate_aa) if val > (r-1)]

for seq in sequences:
    forward_match=tuple(map(CountOccurrences,aa_motifs_forward,seq))
    backward_match=tuple(map(CountOccurrences,aa_motifs_backward))
    total_match=[forward_match[i] if val < (r) else forward_match[i]+backward_match[i] for i, val in enumerate(duplicate_aa)]
    
