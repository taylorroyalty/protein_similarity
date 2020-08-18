import pandas

sequences=pandas.read_table('data/test_sequences.txt')

blosum62=pandas.read_csv('data/blosum62.csv')
blosum62_long=pandas.melt(blosum62,id_vars="aa",var_name="aa2",value_name="score")
seq_align_score=pandas.concat([blosum62_long["aa"]+blosum62_long["aa2"],blosum62_long["score"]],axis=1)


for i, seq1 in sequences.iterrows():
    for j, seq2 in sequences.iloc[i+1:].iterrows():
        seq1_list=list(seq1['sequences'])
        seq2_list=list(seq2['sequences'])
        seq_align_list=tuple(map(''.join,zip(seq1_list,seq2_list)))
        seq_align_obs=pandas.DataFrame(seq_align,columns=["alignment"])
        seq_align_score

