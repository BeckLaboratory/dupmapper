#!/usr/bin/env python3

"""

"""

import argparse
import os
import pandas as pd
import pybedtools
import sys

print('filter_blast.py starting', file=sys.stderr, flush=True)

BLAST_COLS = ['#CHROM', 'POS', 'END', 'TIG', 'TIG_POS', 'TIG_END', 'STRAND', 'EVALUE', 'MISMATCH', 'GAPS', 'PINDENT']

if __name__ == '__main__':

    ap = argparse.ArgumentParser(
        prog='filter_blast.py',
        description='Filter formatted BLAST output against a BED file'
    )

    ap.add_argument('-b', '--bed', dest='bed_filter', required=True, help='BED file to filter out.')
    ap.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true', help='Output progress.')
    ap.add_argument('-t', '--tempdir', dest='temp_dir', default='/tmp', help='Temporary directory for pybedtools')
    ap.add_argument('-f', '--frac', dest='frac_overlap', default=0.5, type=float, help='Remove BLAST records with this proportion (0.0 - 1.0) intersecting filtered loci.')
    ap.add_argument('-c', '--chunksize', dest='chunk_size', default=5000, type=int, help='Process BLAST hits in chunks of this many records. High values use more memory, very low values will incur a performance penalty.')

    args = ap.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)
    pybedtools.set_tempdir(args.temp_dir)

    # Read BED
    if args.verbose:
        print(f'Reading filter: {args.bed_filter}', file=sys.stderr, flush=True)

    bed_filter = pybedtools.BedTool.from_dataframe(
        pd.read_csv(args.bed_filter, sep='\t', header=None, low_memory=False)
    )

    if args.verbose:
        print(f'Done reading filter.', file=sys.stderr,  flush=True)

    # Read dataframe from stdin in chunks
    df_iter = pd.read_csv(
        sys.stdin, iterator=True, chunksize=args.chunk_size,
        sep='\t',
        names=BLAST_COLS,
        header=None
    )

    print_header = True

    chunk_count = 0

    for df in df_iter:

        chunk_count += 1

        if args.verbose:
            print(f'Chunk {chunk_count}: rows = {df.shape[0]}', file=sys.stderr,  flush=True)

        if print_header:
            df.head(0).to_csv(sys.stdout, sep='\t', index=False, header=True)
            print_header = False

        df['POS'] -= 1

        for index, row in df.iterrows():
            if row['POS'] > row['END']:
                (df.loc[index, 'POS'], df.loc[index, 'END']) = (row['END'], row['POS'])

        df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        df.insert(3, 'INDEX', df.index)

        if df.shape[0] == 0:
            continue

        df_filt = pybedtools.BedTool.from_dataframe(df[['#CHROM', 'POS', 'END', 'INDEX']]).intersect(
            b=bed_filter, wa=True, f=args.frac_overlap, v=True
        ).to_dataframe()

        if df_filt.shape[0] == 0:
            continue

        df = df.loc[list(df_filt['name'])]

        del(df['INDEX'])

        if args.verbose:
            print(f'Found {df.shape[0]} records in chunk {chunk_count}', file=sys.stderr,  flush=True)

        df.to_csv(sys.stdout, sep='\t', index=False, header=False)

print('filter_blast.py done', file=sys.stderr, flush=True)
