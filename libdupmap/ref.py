"""
Reference utilities.
"""

import pandas as pd
import re

import svpoplib


def masked_fasta_to_bed(fa_file_name, soft=True, dist=20):
    """
    Get a DataFrame of masked loci in BED format from a FASTA file. Includes hard-masked loci and soft-masked loci
    by default.

    :param fa_file_name: Masked FASTA file name.
    :param soft: Include soft-masked loci if `True`.
    :param dist: Condense records within this number of bases. Set to 0 to disable condensing records.

    :return: A dataframe of masked loci.
    """

    record_list = list()

    if soft:
        mask_pattern = re.compile('[a-zNn]+')
    else:
        mask_pattern = re.compile('[Nn]+')

    for record in svpoplib.seq.fa_to_record_iter(
        fa_file_name, input_format='fasta'
    ):

        match_iter = re.finditer(mask_pattern, str(record.seq))

        try:
            match = next(match_iter)
        except StopIteration():
            continue  # No masked loci, move to the next sequence record

        cur_pos, cur_end = match.span()
        end_dist = cur_end + dist

        if cur_pos < dist:
            cur_pos = 0

        for match in re.finditer(mask_pattern, str(record.seq)):

            pos, end = match.span()

            if pos > end_dist:  # New record
                record_list.append((record.id, cur_pos, cur_end))

                cur_pos = pos
                cur_end = end
                end_dist = end + dist

            elif end > cur_end:  # Extend current record
                cur_end = end
                end_dist = end + dist

        # Move to end if less than dist
        if len(record.seq) - cur_end < dist:
            cur_end = len(record.seq)

        # Append last record
        record_list.append((record.id, cur_pos, cur_end))

    # Merge records to a table
    return pd.DataFrame(
        record_list, columns=['#CHROM', 'POS', 'END']
    ).astype(
        {'POS': int, 'END': int}
    ).sort_values(['#CHROM', 'POS', 'END'])


def merge_regions(df_list, dist=20):

    # Read all regions into one DataFrame
    def read_df_item(df_item):
        if isinstance(df_item, pd.DataFrame):
            df = df_item.copy()

            df.columns = ['#CHROM', 'POS', 'END']

            df = df.astype(
                {'#CHROM': str, 'POS': int, 'END': int}
            )

            return df

        elif isinstance(df_item, str):

            df_item = df_item.strip()

            if not df_item:
                raise ValueError(f'merge_regions(): df_item contains an empty string')

            try:
                libinsblast.util.check_regular_file(df_item)
            except Exception as ex:
                raise ValueError(f'merge_regions(): df_item contains an invalid filename: {str(ex)}')

            return pd.read_csv(
                df_item, sep='\t',
                header=None,
                usecols=(0, 1, 2), names=('#CHROM', 'POS', 'END'),
                comment='#',
                keep_default_na=False
            ).astype(
                {'#CHROM': str, 'POS': int, 'END': int}
            )

        else:
            raise ValueError(f'merge_regions(): df_item contains an item that is not a DataFrame or str: {str(type(df_item))}')

    df = pd.concat(
        [read_df_item(df_item) for df_item in df_list],
        axis=0
    )

    df.sort_values(['#CHROM', 'POS'], inplace=True)

    # Collapse
    record_list = list()

    chrom = None
    pos = None
    end = None
    end_dist = None

    for index, row in df.iterrows():

        if row['#CHROM'] != chrom:
            # Handle chromosome change

            if chrom is not None:
                record_list.append((chrom, pos, end))

            chrom = row['#CHROM']
            pos = row['POS']
            end = row['END']
            end_dist = end + dist

            continue

        if row['POS'] > end_dist:
            # Current record is out of the last one

            record_list.append((chrom, pos, end))

            pos = row['POS']
            end = row['END']
            end_dist = end + dist

        elif row['END'] > end:
            # Extend the current record

            end = row['END']
            end_dist = end + dist

    # Last record
    if chrom is not None:
        record_list.append((chrom, pos, end))

    # Merge
    return pd.DataFrame(
        record_list, columns=['#CHROM', 'POS', 'END']
    ).astype(
        {'POS': int, 'END': int}
    )
