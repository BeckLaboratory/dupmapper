"""
Filtering routines for BLAST hits.
"""

def filter_ins(df):
    # Filter for BLAST hits between the breakpoints
    df['BLAST_MIN'] = list(blast_min[df['TIG']])
    df['BLAST_MAX'] = list(blast_max[df['TIG']])

    df = df.loc[
        (df['END'] > df['BLAST_MIN']) &
        (df['POS'] < df['BLAST_MAX'])
        ]

    # Save record
    if df.shape[0] > 0:
        del (df['BLAST_MIN'])
        del (df['BLAST_MAX'])
        df_list.append(df)