# BLAST rules

# Write BLAST xlsx
rule eshift_td_anno_insdel_xlsx:
    input:
        bed='results/{sample}/blast/td_sv_{svtype}.bed.gz'
    output:
        xlsx='results/{sample}/blast/td_sv_{svtype}.xlsx'
    run:

        pd.read_csv(
            input.bed, sep='\t'
        ).to_excel(
            output.xlsx, index=False
        )

# Annotate TDs using BLAST hits.
rule eshift_td_anno:
    input:
        bed='temp/{sample}/blast/filtered_td_sv_{svtype}.bed.gz',
        bed_sv='temp/{sample}/input/sv_{svtype}.bed.gz'
    output:
        bed='results/{sample}/blast/td_sv_{svtype}.bed.gz'
    params:
        min_len=30
    threads: 1
    run:

        # Read BLAST records
        df = pd.read_csv(input.bed,sep='\t')

        # Read variants
        df_sv = pd.read_csv(
            input.bed_sv, sep='\t',
            usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'),
            index_col='ID'
        )

        df_sv = df_sv.loc[sorted(set(df['TIG']))]

        # Select BLAST records on the same chromosome as the SV
        if df.shape[0] > 0:
            df = df.loc[df.apply(lambda row: row['#CHROM'] == df_sv.loc[row['TIG'], '#CHROM'], axis=1)].copy()

        # Set BLAST length (reference and tig space)
        df['LEN'] = df['END'] - df['POS']
        df['TIG_LEN'] = df['TIG_END'] - df['TIG_POS']

        df = df.loc[df['LEN'] >= params.min_len]

        # Filter redundant hits in reference space
        df.sort_values('LEN', ascending=False, inplace=True)

        index_set = set()

        for var_id, subdf in df.groupby('TIG'):

            if subdf.shape[0] == 1:
                index_set.add(subdf.iloc[0].name)
                continue

            itree = intervaltree.IntervalTree()

            for index, row in subdf.iterrows():
                keep = True

                interval_set = itree[row['POS']:row['END']]

                for interval in interval_set:
                    if (
                        min([row['END'], interval.end]) - max([row['POS'], interval.begin])
                    ) / row['LEN'] > 0.8:
                        keep = False

                if keep:
                    itree[row['POS']:row['END']] = index
                    index_set.add(index)

        df = df.loc[sorted(index_set)]

        if df.shape[0] > 0:
            # Filter redundant hits in contig (SV sequence) space
            df.sort_values('TIG_LEN', ascending=False, inplace=True)

            index_set = set()

            for var_id, subdf in df.groupby('TIG'):

                if subdf.shape[0] == 1:
                    index_set.add(subdf.iloc[0].name)
                    continue

                itree = intervaltree.IntervalTree()

                for index, row in subdf.iterrows():
                    keep = True

                    interval_set = itree[row['TIG_POS']:row['TIG_END']]

                    for interval in interval_set:
                        if (
                            min([row['TIG_END'], interval.end]) - max([row['TIG_POS'], interval.begin])
                        ) / row['TIG_LEN'] > 0.8:
                            keep = False

                    if keep:
                        itree[row['TIG_POS']:row['TIG_END']] = index
                        index_set.add(index)

            df = df.loc[sorted(index_set)]

            # Calculate span (90% of SVLEN aligned by BLAST, less than 10% difference in overall span vs SVLEN)
            df_sv = df_sv.loc[sorted(set(df['TIG']))]

            blast_bp = df.groupby('TIG').apply(libdupmap.blast.td_count_span_bp)

            if wildcards.svtype == 'ins':
                blast_span = df.groupby('TIG').apply(lambda subdf:
                    np.max(subdf['END']) - np.min(subdf['POS'])
                )

            else:

                blast_span = pd.Series(0, index=set(df['TIG']))

                for var_id, subdf in df.groupby('TIG'):
                    pos = df_sv.loc[var_id, 'POS']
                    end = df_sv.loc[var_id, 'END']

                    subdf = subdf.copy()

                    subdf['SIDE'] = 'N'  # Unknown side of the deletion (filled in below)

                    for index in subdf.index:
                        b_pos = subdf.loc[index, 'POS']
                        b_end = subdf.loc[index, 'END']

                        if b_pos < pos and b_end < end:
                            # Trim record that extends into the deletion from the left
                            subdf.loc[index, 'END'] = min([b_end, pos])
                            subdf.loc[index, 'SIDE'] = 'L'

                        elif b_pos > pos and b_end > end:
                            # Trim record that extends into the deletion from right left
                            subdf.loc[index, 'POS'] = max([b_pos, end])
                            subdf.loc[index, 'SIDE'] = 'R'

                        else:
                            # Remove BLAST records if they span the deletion
                            subdf = subdf.loc[[val != index for val in subdf.index]]

                    if subdf.shape[0] == 0:
                        continue

                    # Remove records that were collapsed to nothing
                    subdf = subdf.loc[subdf['END'] - subdf['POS'] > 0]

                    if subdf.shape[0] == 0:
                        continue

                    # Recalculate SPAN outside DEL
                    blast_span[var_id] = np.sum(
                        subdf.groupby('SIDE').apply(lambda subsubdf:
                            np.max(subsubdf['END']) - np.min(subsubdf['POS'])
                        )
                    )

                blast_span = blast_span.loc[blast_span > 0]

            prop_span = np.abs(blast_span - df_sv['SVLEN']) / df_sv['SVLEN']

            prop_span = prop_span.loc[(prop_span < 0.1) & (blast_bp / df_sv['SVLEN'] > 0.9)]

            # Subset dataframes to actual TDs
            id_set = set(prop_span.index)

            df = df.loc[df['TIG'].apply(lambda val: val in id_set)].copy()
            df_sv = df_sv.loc[sorted(id_set)]

            # Determine if records are left or right of insertion site
            if wildcards.svtype == 'ins':
                df['BREAK_SIDE'] = df.apply(lambda row:
                    'L' if ((row['POS'] + row['END']) / 2) < df_sv.loc[row['TIG'], 'POS'] else 'R',
                    axis=1
                )
            else:
                df['BREAK_SIDE'] = df.apply(lambda row:
                    'L' if ((row['POS'] + row['END']) / 2) < ((df_sv.loc[row['TIG'], 'POS'] + df_sv.loc[row['TIG'], 'POS']) / 2) else 'R',
                    axis=1
                )

            # make final table
            df_list = list()

            for var_id, subdf in df.groupby('TIG'):

                row_dtab = df_sv.loc[var_id]

                n_l = np.sum(subdf['BREAK_SIDE'] == 'L')
                n_r = np.sum(subdf['BREAK_SIDE'] == 'R')

                if n_l == 1:
                    row = subdf.loc[subdf['BREAK_SIDE'] == 'L'].squeeze()
                    prop_l = 1.0 if row['END'] < row_dtab['POS'] else (row['END'] - row_dtab['POS']) / row['LEN']
                    sv_l = f'{row["TIG_POS"]}-{row["TIG_END"]}'
                else:
                    prop_l = np.nan
                    sv_l = np.nan

                if n_r == 1:
                    row = subdf.loc[subdf['BREAK_SIDE'] == 'R'].squeeze()
                    prop_r = 1.0 if row['POS'] > row_dtab['POS'] else (row_dtab['POS'] - row['POS']) / row['LEN']
                    sv_r = f'{row["TIG_POS"]}-{row["TIG_END"]}'
                else:
                    prop_r = np.nan
                    sv_r = np.nan

                map_pos = np.min(subdf['POS'])
                map_end = np.max(subdf['END'])

                tig_pos = np.min(subdf['TIG_POS'])
                tig_end = np.min(subdf['TIG_POS'])

                df_list.append(
                    pd.Series(
                        [
                            var_id,
                            map_pos, map_end, map_end - map_pos,
                            tig_pos, tig_end, tig_end - tig_pos,
                            n_l, n_r, prop_l, prop_r,
                            sv_l, sv_r
                        ],
                        index=[
                            'ID',
                            'MAP_POS', 'MAP_END', 'MAP_LEN',
                            'SV_MAP_POS', 'SV_MAP_END', 'SV_MAP_LEN',
                            'MAP_L_N', 'MAP_R_N', 'MAP_L_PROP', 'MAP_R_PROP',
                            'SV_L', 'SV_R'
                        ]
                    )
                )

            df_td = pd.concat(
                [
                    df_sv,
                    pd.concat(df_list, axis=1).T.set_index('ID')
                ], axis=1
            )

        else:
            df_td = pd.DataFrame(
                [],
                columns=[
                            'ID',
                            'MAP_POS', 'MAP_END', 'MAP_LEN',
                            'SV_MAP_POS', 'SV_MAP_END', 'SV_MAP_LEN',
                            'MAP_L_N', 'MAP_R_N', 'MAP_L_PROP', 'MAP_R_PROP',
                            'SV_L', 'SV_R'
                        ]
            )

        # Write
        df_td.to_csv(output.bed, sep='\t', index=True, compression='gzip')

# Apply basic filters - BLAST hits must land proximally to SV breakpoints
rule dmap_td_blast_prox:
    input:
        bed='results/{sample}/blast/blast_full.bed.gz',
        bed_sv='temp/{sample}/input/sv_{svtype}.bed.gz'
    output:
        bed=temp('temp/{sample}/blast/filtered_td_sv_{svtype}.bed.gz')
    params:
        chunksize=50000
    wildcard_constraints:
        svtype='ins|del'
    threads: 1
    run:

        # Read variant BED
        df_sv = pd.read_csv(input.bed_sv, sep='\t')
        df_sv = df_sv.loc[(df_sv['SVTYPE'] == wildcards.svtype.upper()) & (df_sv['PART_BLAST'] >= 0)]

        id_set = df_sv['ID']

        # Set range where BLAST hits must be found
        df_sv['PAD'] = df_sv['SVLEN'].apply(lambda val: max([val * config['td_pad'], config['td_pad_min']])).astype(int)
        df_sv['BLAST_MIN'] = df_sv['POS'] - df_sv['PAD']
        df_sv['BLAST_MAX'] = df_sv['END'] + df_sv['PAD']

        blast_min = df_sv.set_index('ID')['BLAST_MIN']
        blast_max = df_sv.set_index('ID')['BLAST_MAX']

        sv_pos = df_sv.set_index('ID')['POS']
        sv_end = df_sv.set_index('ID')['END']

        # Filter BLAST hits
        df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=params.chunksize)

        df_list = list()

        header_list = None

        for df in df_iter:

            if header_list is None:
                header_list = list(df.columns)

            df = df.loc[df['TIG'].isin(id_set)]

            if df.shape[0] == 0:
                continue

            # Filter for BLAST hits between the breakpoints
            df['BLAST_MIN'] = list(blast_min[df['TIG']])
            df['BLAST_MAX'] = list(blast_max[df['TIG']])

            df = df.loc[
                (df['END'] > df['BLAST_MIN']) &
                (df['POS'] < df['BLAST_MAX'])
            ]

            # Remove direct deletion region hits
            if wildcards.svtype == 'del':

                df['SV_POS'] = list(sv_pos[df['TIG']])
                df['SV_END'] = list(sv_end[df['TIG']])

                ro = df.apply(lambda row:
                svpoplib.variant.reciprocal_overlap(row['POS'],row['END'],row['SV_POS'],row['SV_END']),
                    axis=1
                )

                df = df.loc[ro < 0.5]

                del(df['SV_POS'])
                del(df['SV_END'])

            # Save record
            if df.shape[0] > 0:
                del(df['BLAST_MIN'])
                del(df['BLAST_MAX'])
                df_list.append(df)

        if df_list:
            df = pd.concat(df_list, axis=0)
        else:
            df = pd.DataFrame([], columns=header_list)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# Merge full BLAST output.
rule dmap_blast_merge_full:
    input:
        bed=expand('temp/{{sample}}/blast/part/part_{part}.bed.gz', part=range(config['partitions']))
    output:
        bed='results/{sample}/blast/blast_full.bed.gz'
    threads: 1
    run:

        write_header = True

        with gzip.open(output.bed, 'wt') as out_file:
            for in_file_name in input.bed:
                pd.read_csv(
                    in_file_name, sep='\t'
                ).to_csv(
                    out_file, sep='\t', index=False, header=write_header
                )

                write_header = False

# Run BLAST.
rule dmap_blast_run:
    input:
        fa='temp/{sample}/input/fa_part/part_{part}.fa',
        bed_filter='data/ref/filter.bed.gz'
    output:
        bed=temp('temp/{sample}/blast/part/part_{part}.bed.gz')
    threads: 32
    params:
        temp_dir=lambda wildcards: os.path.join(config['temp_dir'], 'blast_{wildcards.part}'),
        blast_db=config['blast_db'],
        blast_param=config['blast_param'],
        chunksize=5000,
        filter_blast=os.path.join(DMAP_DIR, 'scripts', 'filter_blast.py')
    shell:
        """blastn -db {params.blast_db} """
            """-num_threads {threads} """
            """-outfmt "6 saccver sstart send qseqid qstart qend sstrand evalue mismatch gaps pident" """
            """{params.blast_param} """
            """-query {input.fa} | """
        """{params.filter_blast} -b {input.bed_filter} -t {params.temp_dir} -c {params.chunksize} -v | """
        """gzip > {output.bed}"""
