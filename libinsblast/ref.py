"""
Reference utilities.
"""

import pandas as pd
import re

import svpoplib


def masked_fasta_to_bed(fa_file_name, soft=True):

    record_list = list()

    if soft:
        mask_pattern = re.compile('[a-zNn]+')
    else:
        mask_pattern = re.compile('[Nn]+')

    for record in svpoplib.seq.fa_to_record_iter(
        fa_file_name, input_format='fasta'
    ):

        # Traverse reference with a state machine:
        # * U: Unique (not masked)
        # * M: Masked

        for match in re.finditer(mask_pattern, str(seq.seq)):
            pos, end = match.span()
            record_list.append((record.id, pos, end))

    return pd.DataFrame(
        record_list, columns=['#CHROM', 'POS', 'END']
    ).astype(
        {'POS': int, 'END': int}
    )
