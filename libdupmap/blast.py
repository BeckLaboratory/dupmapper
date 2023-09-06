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


def td_count_span_bp(subdf):
    """
    Count the number of bases covered by BLAST alignments.
    """

    subdf = subdf.sort_values(['POS', 'END'])

    last_pos = subdf.iloc[0]['POS']
    last_end = subdf.iloc[0]['END']

    if subdf.shape[0] == 1:
        return last_end - last_pos

    span_bp = 0

    for index, row in subdf.iloc[1:].iterrows():
        if row['POS'] > last_end:
            span_bp += last_end - last_pos
            last_pos = row['POS']

        last_end = row['END']

    span_bp += last_end - last_pos

    return span_bp