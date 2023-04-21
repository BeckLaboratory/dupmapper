"""
Map insertion sequences back to the reference with BLAST to find patterns of duplications.
"""

import os
import sys

sys.path.append(os.path.join('dep', 'svpop'))
sys.path.append(os.path.join('dep', 'svpop', 'dep'))
sys.path.append(os.path.join('dep', 'svpop', 'dep', 'ply'))

import libinsblast
import svpoplib


#
# Definitions
#


#
# Init
#

# Read and check configuration file

config_file_name = libinsblast.pipeline.find_config(config)

if config_file_name is not None:
    configfile: config_file_name

libinsblast.pipeline.check_config(config)


#
# Rules
#

localrules: eshift_td_remap_all

rule eshift_td_remap_all:
    input:
        stats=expand('results/td/stats/{subset}/td_stats_sv_insdel.tsv.gz', subset=('notrsd', 'sdnotr')),
        boxplot_max=expand('results/td/fig/{subset}/box_max_offset/notrsd/td_box_max-bp_sv_insdel.png', subset=('notrsd', 'sdnotr'))



eshift_td_blast_mem = {  # Increase memory for batches that need it
    ('notrsd', 'ins', '13'): '--mem=131072 -c 16 -t 48:00'
}

# eshift_merge_td_boxplot
#
# Make a boxplot of the breakpoint offsets for TDs.
rule eshift_merge_td_boxplot:
    input:
        tsv_ins='results/key/dtab/dtab_sv_ins.tsv.gz',
        tsv_del='results/key/dtab/dtab_sv_del.tsv.gz',
        bed_td='results/td/anno/{subset}/td_anno_sv_insdel.bed.gz',
        tsv_stat='results/td/stats/{subset}/td_stats_sv_insdel.tsv.gz'
    output:
        png='results/td/fig/{subset}/box_max_offset/notrsd/td_box_max-bp_sv_insdel.png',
        pdf='results/td/fig/{subset}/box_max_offset/notrsd/td_box_max-bp_sv_insdel.pdf',
        svg='results/td/fig/{subset}/box_max_offset/notrsd/td_box_max-bp_sv_insdel.svg',
        png_low='results/td/fig/{subset}/box_max_offset/notrsd/td_box_max-bp_sv_insdel_lowres.png'
    params:
        cluster_opts='--mem=2048 -t 00:30'
    run:

        subset = 'notrsd'

        FIG_SIZE = (7, 7)
        FIG_DPI = 300

        BOX_COLOR = {
            'ntd': 'darkgreen',
            'td': 'mediumorchid'
        }

        FP_STU_MIN = 0.01

        # Read breakpoint distribution
        df_ins = pd.read_csv(input.tsv_ins, sep='\t')
        df_del = pd.read_csv(input.tsv_del, sep='\t')

        df = pd.concat([df_ins, df_del], axis=0)

        df = subset_dtab(df, wildcards.subset)

        del(df_ins)
        del(df_del)

        df = df.loc[(df['AC'] > 1) & (df['OFFSET_MAX'] > 0)].copy()

        td_id_set = set(
            pd.read_csv(input.bed_td, sep='\t', usecols=('ID', ))['ID']
        )

        df['TD'] = df['ID'].apply(lambda val: val in td_id_set)

        # Read stats for significance
        df_stat = pd.read_csv(input.tsv_stat, sep='\t', index_col=('SUBSET', 'SVTYPE'))
        df_stat = df_stat.loc[df_stat['SHIFT_ONLY'] & (df_stat['VARTYPE'] == 'sv')]

        # Adjust for deletions
        max_offset = np.max(df.loc[df['SVTYPE'] == 'INS', 'OFFSET_MAX'])
        del_chopped_ntd = np.sum((df['SVTYPE'] == 'DEL') & (df['OFFSET_MAX'] > max_offset) & ~ df['TD'])
        del_chopped_td = np.sum((df['SVTYPE'] == 'DEL') & (df['OFFSET_MAX'] > max_offset) & df['TD'])

        df = df.loc[(df['SVTYPE'] == 'INS') | (df['OFFSET_MAX'] <= max_offset)]

        # Make figure
        ax_dict = dict()
        fig, (ax_dict['ins'], ax_dict['del']) = plt.subplots(1, 2, figsize=FIG_SIZE, dpi=FIG_DPI, sharey='all')

        for svtype in ('ins', 'del'):

            ax = ax_dict[svtype]

            df_svtype = df.loc[df['SVTYPE'] == svtype.upper()]

            bplot = ax.boxplot(
                [
                    df_svtype.loc[~ df_svtype['TD'], 'OFFSET_MAX'],
                    df_svtype.loc[df_svtype['TD'], 'OFFSET_MAX']
                ],
                sym='ko',
                widths=(0.75, 0.75),
                notch=True,
                patch_artist=True
            )

            ax.set_title(
                SVTYPE_NAME_PLU[svtype],
                fontfamily='arial', size=18
            )

            ax.set_xticks(
                [1, 2],
                ['non-TD', 'TD'],
                fontfamily='arial',
                size=16
            )

            # Box aestetics
            bplot['boxes'][0].set_facecolor(BOX_COLOR['ntd'])
            bplot['boxes'][1].set_facecolor(BOX_COLOR['td'])

            bplot['medians'][0].set_color('black')
            bplot['medians'][1].set_color('black')

            bplot['medians'][0].set_linewidth(2)
            bplot['medians'][1].set_linewidth(2)

            # Alpha fliers
            bplot['fliers'][0].set_alpha(0.6)
            bplot['fliers'][1].set_alpha(0.6)

            bplot['fliers'][0].set_color(BOX_COLOR['ntd'])
            bplot['fliers'][1].set_color(BOX_COLOR['td'])

            # Spines and axes
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)

            if svtype != 'ins':
                ax.tick_params(axis='y', which='both', left=False)

        # Points off axis
        y_max = ax_dict['ins'].get_ylim()[1]

        arrow_height = y_max * 0.03
        arrow_width = 0.1
        arrow_color = 'firebrick'

        if del_chopped_td > 0 or del_chopped_ntd > 0:
            ax = ax_dict['del']

            if del_chopped_ntd > 0:
                ax.arrow(
                    1, y_max, 0, arrow_height, fill=True,
                    head_length=arrow_height, head_width=arrow_width,
                    length_includes_head=True,
                    color=arrow_color
                )

                f'{del_chopped_ntd:,d}'

                ax.text(
                    1, y_max + arrow_height + y_max * 0.01,
                    f'{del_chopped_ntd:,d}',
                    size=12,
                    color=arrow_color,
                    ha='center',
                    fontfamily='arial'
                )

            if del_chopped_td > 0:
                ax.arrow(
                    2, y_max, 0, arrow_height, fill=True,
                    head_length=arrow_height, head_width=arrow_width,
                    length_includes_head=True,
                    color=arrow_color
                )

                f'{del_chopped_ntd:,d}'

                ax.text(
                    2, y_max + arrow_height + y_max * 0.01,
                    f'{del_chopped_td:,d}',
                    size=12,
                    color=arrow_color,
                    ha='center',
                    fontfamily='arial'
                )

            # Adjust space for significance bars
            y_max *= 1.07

        # Significance bars
        sig_height = y_max * 1.05
        sig_leg_len = y_max * 0.02

        for svtype in ('ins', 'del'):
            ax = ax_dict[svtype]

            stat_row = df_stat.loc[(subset, svtype)].squeeze()

            sig = stat_row['T_P']

            if sig > 0.01:
                sig_txt = 'n.s.'
            elif sig > 0.001:
                sig_txt = '*'
            elif sig > 0.0001:
                sig_txt = '**'
            else:
                sig_txt = '***'

            ax.plot(
                [1, 1, 2, 2],
                [sig_height - sig_leg_len, sig_height, sig_height, sig_height - sig_leg_len],
                linewidth=0.78, color='black'
            )

            if sig_txt != 'n.s.':
                ax.text(1.5, sig_height - sig_leg_len, sig_txt, ha='center', va='bottom', c='black', fontweight='bold', fontfamily='arial', size=20)
            else:
                ax.text(1.5, sig_height, sig_txt, ha='center', va='bottom', c='black', fontweight='bold', fontfamily='arial', size=16)

        # Write
        fig.savefig(output.png)
        fig.savefig(output.pdf)
        fig.savefig(output.svg)
        fig.savefig(output.png_low, dpi=76)

# eshift_td_stats
#
# Get shifted stats for TDs.
rule eshift_td_stats:
    input:
        bed_td_ins='results/td/anno/{subset}/td_anno_sv_ins.bed.gz',
        bed_td_del='results/td/anno/{subset}/td_anno_sv_del.bed.gz',
        tsv_ins='results/key/dtab/dtab_sv_ins.tsv.gz',
        tsv_del='results/key/dtab/dtab_sv_del.tsv.gz'
    output:
        tsv='results/td/stats/{subset}/td_stats_sv_insdel.tsv.gz',
        xlsx='results/td/stats/{subset}/td_stats_sv_insdel.xlsx'
    params:
        cluster_opts='--mem=1024 -t 00:15'
    run:

        vartype = 'sv'

        # Get stats
        df_stat_list = list()

        for svtype in SVTYPES:

            df = pd.read_csv(input[f'tsv_{svtype}'], sep='\t')
            df = subset_dtab(df, wildcards.subset)
            dtab_id_set = set(df['ID'])

            df_td = pd.read_csv(input[f'bed_td_{svtype}'], sep='\t')
            df_td = df_td.loc[df_td['ID'].apply(lambda val: val in dtab_id_set)]
            td_id_set = set(df_td['ID'])

            df_td['ALN_ROT'] = ~ (pd.isnull(df_td['SV_L']) | pd.isnull(df_td['SV_L']))

            df['TD'] = df['ID'].apply(lambda val: val in td_id_set)

            df_ac1 = df.loc[df['AC'] > 1]

            max_td = df_ac1.loc[df_ac1['TD'], 'OFFSET_MAX']
            max_ntd = df_ac1.loc[~df_ac1['TD'], 'OFFSET_MAX']

            f_val, f_p = f_test(max_td, max_ntd)

            t_stu = f_p >= F_P_SIG

            t, t_p = scipy.stats.ttest_ind(
                max_td,
                max_ntd,
                equal_var=t_stu,
                alternative='two-sided'
            )

            cohen_d = get_cohen_d(max_td, max_ntd)

            fisher_or, fisher_p = scipy.stats.fisher_exact(
                [
                    [np.sum(max_td == 0), np.sum(max_td > 0)],
                    [np.sum(max_ntd == 0), np.sum(max_ntd > 0)]
                ]
            )

            df_stat_list.append(
                pd.Series(
                    [
                        wildcards.subset, vartype, svtype,
                        df_td.shape[0], len(dtab_id_set - td_id_set), df_td.shape[0] / len(dtab_id_set),
                        np.sum(df_td['ALN_ROT']), np.sum(df_td['ALN_ROT']) / df_td.shape[0],
                        max_td.shape[0], max_ntd.shape[0], max_td.shape[0] / (max_td.shape[0] + max_ntd.shape[0]),
                        np.sum(max_td > 0), np.sum(max_td > 0) / max_td.shape[0],
                        np.mean(max_td), np.mean(max_ntd),
                        np.std(max_td), np.std(max_ntd),
                        f_val, f_p,
                        t, t_p, 'STU' if t_stu else 'WEL',
                        cohen_d,
                        fisher_or, fisher_p
                    ],
                    index=[
                        'SUBSET', 'VARTYPE', 'SVTYPE',
                        'N_TD', 'N_NTD', 'PROP_TD',
                        'ALN_ROT', 'ALN_ROT_PROP',
                        'N_TD_AC1', 'N_NTD_AC1', 'PROP_TD_AC1',
                        'N_TD_SHIFT', 'PROP_TD_SHIFT',
                        'MEAN_TD_MAX', 'MEAN_NTD_MAX',
                        'SD_TD_MAX', 'SD_NTD_MAX',
                        'F_MAX', 'F_P_MAX',
                        'T_MAX', 'T_P_MAX', 'T_TEST_MAX',
                        'COHEN_D',
                        'FET_OR', 'FET_P'
                    ]
                )
            )

        # Merge and write
        df_stat = pd.concat(df_stat_list, axis=1).T

        for col in ('N_TD', 'N_NTD'):
            df_stat[col] = df_stat[col].astype(int)

        df_stat.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
        df_stat.to_excel(output.xlsx, index=False)

# eshift_td_anno_insdel
#
# Merge INS/DEL
rule eshift_td_anno_insdel:
    input:
        bed_ins='results/td/anno/{subset}/td_anno_sv_ins.bed.gz',
        bed_del='results/td/anno/{subset}/td_anno_sv_del.bed.gz'
    output:
        bed='results/td/anno/{subset}/td_anno_sv_insdel.bed.gz',
        xlsx='results/td/anno/{subset}/td_anno_sv_insdel.xlsx'
    params:
        cluster_opts='--mem=512 -t 00:10'
    run:

        df = pd.concat(
            [
                pd.read_csv(input.bed_ins, sep='\t'),
                pd.read_csv(input.bed_del, sep='\t')
            ], axis=0
        )

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
        df.to_excel(output.xlsx, index=False)

# eshift_td_anno_ins
#
# Annotate INS TDs using BLAST hits.
rule eshift_td_anno_ins:
    input:
        bed='results/td/blast/filtered/{subset}/filtered_td_blast_sv_{svtype}.bed.gz',
        tsv_dtab='results/key/dtab/dtab_sv_{svtype}.tsv.gz'
    output:
        bed='results/td/anno/{subset}/td_anno_sv_{svtype}.bed.gz',
        xlsx='results/td/anno/{subset}/td_anno_sv_{svtype}.xlsx'
    params:
        cluster_opts='--mem=1024 -t 00:30'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        df_dtab = pd.read_csv(
            input.tsv_dtab, sep='\t',
            usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'),
            index_col='ID'
        )

        df_dtab = df_dtab.loc[sorted(set(df['TIG']))]

        # Filter df
        df = df.loc[df.apply(lambda row: row['#CHROM'] == df_dtab.loc[row['TIG'], '#CHROM'], axis=1)].copy()

        df['LEN'] = df['END'] - df['POS']
        df['TIG_LEN'] = df['TIG_END'] - df['TIG_POS']

        df = df.loc[df['LEN'] >= 30]

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
        df_dtab = df_dtab.loc[sorted(set(df['TIG']))]

        blast_bp = df.groupby('TIG').apply(eshift_td_count_span_bp)

        if wildcards.svtype == 'ins':
            blast_span = df.groupby('TIG').apply(lambda subdf:
                np.max(subdf['END']) - np.min(subdf['POS'])
            )

        else:

            blast_span = pd.Series(0, index=set(df['TIG']))

            for var_id, subdf in df.groupby('TIG'):
                pos = df_dtab.loc[var_id, 'POS']
                end = df_dtab.loc[var_id, 'END']

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

        prop_span = np.abs(blast_span - df_dtab['SVLEN']) / df_dtab['SVLEN']

        prop_span = prop_span.loc[(prop_span < 0.1) & (blast_bp / df_dtab['SVLEN'] > 0.9)]

        # Subset dataframes to actual TDs
        id_set = set(prop_span.index)

        df = df.loc[df['TIG'].apply(lambda val: val in id_set)].copy()
        df_dtab = df_dtab.loc[sorted(id_set)]

        # Determine if records are left or right of insertion site
        if wildcards.svtype == 'ins':
            df['BREAK_SIDE'] = df.apply(lambda row:
                'L' if ((row['POS'] + row['END']) / 2) < df_dtab.loc[row['TIG'], 'POS'] else 'R',
                axis=1
            )
        else:
            df['BREAK_SIDE'] = df.apply(lambda row:
                'L' if ((row['POS'] + row['END']) / 2) < ((df_dtab.loc[row['TIG'], 'POS'] + df_dtab.loc[row['TIG'], 'POS']) / 2) else 'R',
                axis=1
            )

        # make final table
        df_list = list()

        for var_id, subdf in df.groupby('TIG'):

            row_dtab = df_dtab.loc[var_id]

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
                df_dtab,
                pd.concat(df_list, axis=1).T.set_index('ID')
            ], axis=1
        )

        # Write
        df_td.to_csv(output.bed, sep='\t', index=True, compression='gzip')
        df_td.to_excel(output.xlsx, index=True)


# eshift_td_blast_filter_basic_del
#
# Apply basic filters: DEL
rule eshift_td_blast_filter_basic_del:
    input:
        bed='results/td/blast/full/{subset}/full_td_blast_sv_del.bed.gz',
        tsv_dtab='results/key/dtab/dtab_sv_del.tsv.gz'
    output:
        bed='results/td/blast/filtered/{subset}/filtered_td_blast_sv_del.bed.gz'
    params:
        cluster_opts='--mem=2048 -t 48:00'
    run:

        svtype = 'del'

        # Read merged variants
        df_dtab = subset_dtab(
            pd.read_csv(input.tsv_dtab, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'IS_TR', 'IS_SD')),
            wildcards.subset, svtype
        )

        # Set range where BLAST hits must be found
        df_dtab['PAD'] = df_dtab['SVLEN'].apply(lambda val: max([val * .1, 100])).astype(int)
        df_dtab['BLAST_MIN'] = df_dtab['POS'] - df_dtab['PAD']
        df_dtab['BLAST_MAX'] = df_dtab['END'] + df_dtab['PAD']

        blast_min = df_dtab.set_index('ID')['BLAST_MIN']
        blast_max = df_dtab.set_index('ID')['BLAST_MAX']
        sv_pos = df_dtab.set_index('ID')['POS']
        sv_end = df_dtab.set_index('ID')['END']

        # Filter BLAST hits
        df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=50000)

        df_list = list()

        for df in df_iter:

            # Filter for BLAST hits between the breakpoints
            df['BLAST_MIN'] = list(blast_min[df['TIG']])
            df['BLAST_MAX'] = list(blast_max[df['TIG']])

            df = df.loc[
                (df['END'] > df['BLAST_MIN']) &
                (df['POS'] < df['BLAST_MAX'])
            ]

            # Remove direct deletion region hits
            df['SV_POS'] = list(sv_pos[df['TIG']])
            df['SV_END'] = list(sv_end[df['TIG']])

            ro = df.apply(lambda row:
                svpoplib.variant.reciprocal_overlap(row['POS'], row['END'], row['SV_POS'], row['SV_END']),
                axis=1
            )

            df = df.loc[ro < 0.5]

            # Save
            if df.shape[0] > 0:
                del(df['BLAST_MIN'])
                del(df['BLAST_MAX'])
                del(df['SV_POS'])
                del(df['SV_END'])
                df_list.append(df)

        df = pd.concat(df_list, axis=0)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# eshift_td_blast_filter_basic_ins
#
# Apply basic filters: INS
rule eshift_td_blast_filter_basic_ins:
    input:
        bed='results/td/blast/full/{subset}/full_td_blast_sv_ins.bed.gz',
        tsv_dtab='results/key/dtab/dtab_sv_ins.tsv.gz'
    output:
        bed='results/td/blast/filtered/{subset}/filtered_td_blast_sv_ins.bed.gz'
    params:
        cluster_opts='--mem=2048 -t 48:00'
    run:

        svtype = 'ins'

        # Read merged variants
        df_dtab = subset_dtab(
            pd.read_csv(input.tsv_dtab, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'IS_TR', 'IS_SD')),
            wildcards.subset, svtype
        )

        # Set range where BLAST hits must be found
        df_dtab['PAD'] = df_dtab['SVLEN'].apply(lambda val: max([val * .1, 100])).astype(int)
        df_dtab['BLAST_MIN'] = df_dtab['POS'] - df_dtab['PAD']
        df_dtab['BLAST_MAX'] = df_dtab['END'] + df_dtab['PAD']

        blast_min = df_dtab.set_index('ID')['BLAST_MIN']
        blast_max = df_dtab.set_index('ID')['BLAST_MAX']

        # Filter BLAST hits
        df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=50000)

        df_list = list()

        for df in df_iter:

            # Filter for BLAST hits between the breakpoints
            df['BLAST_MIN'] = list(blast_min[df['TIG']])
            df['BLAST_MAX'] = list(blast_max[df['TIG']])

            df = df.loc[
                (df['END'] > df['BLAST_MIN']) &
                (df['POS'] < df['BLAST_MAX'])
            ]

            # Save record
            if df.shape[0] > 0:
                del(df['BLAST_MIN'])
                del(df['BLAST_MAX'])
                df_list.append(df)

        df = pd.concat(df_list, axis=0)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# eshift_td_blast_merge_full
#
# Merge full BLAST output.
rule eshift_td_blast_merge_full:
    input:
        bed=expand('temp/td/blast/raw_blast/batch/{{subset}}/sv_{{svtype}}/batch_{batch}.bed.gz', batch=range(TD_FA_BATCH))
    output:
        bed='results/td/blast/full/{subset}/full_td_blast_sv_{svtype}.bed.gz'
    params:
        cluster_opts='--mem=2048 -t 08:00'
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

# eshift_td_blast
#
# Run BLAST.
rule eshift_td_blast:
    input:
        fa='temp/td/remap/fa_batch/{subset}/sv_{svtype}/batch_{batch}.fa',
        fa_filter='results/td/blast/filter/rmsk_tr_filter.bed.gz'
    output:
        bed=temp('temp/td/blast/raw_blast/batch/{subset}/sv_{svtype}/batch_{batch}.bed.gz')
    params:
        threads=16,
        cluster_opts=lambda wildcards: eshift_td_blast_mem.get((wildcards.subset, wildcards.svtype, wildcards.batch), '--mem=32768 -c 16 -t 24:00'),
        temp_dir='/fastscratch/audanp/merge_strategies/ebert_shifted/temp/td_blast_pybed/sv_{svtype}_{batch}'
    shell:
        """{blastn} -db {BLAST_DB} """
            """-num_threads {params.threads} """
            """-outfmt "6 saccver sstart send qseqid qstart qend sstrand evalue mismatch gaps pident" """
            """-word_size 16 """
            """-perc_identity 95 """
            """-query {input.fa} | """
        """scripts/filter_blast.py -b {input.fa_filter} -t {params.temp_dir} -v | """
        """gzip > {output.bed}"""


# eshift_td_blast_input_fa
#
# Split FASTA
rule eshift_td_blast_input_fa:
    input:
        tsv_dtab='results/key/dtab/dtab_sv_{svtype}.tsv.gz'
    output:
        fa=expand('temp/td/remap/fa_batch/{{subset}}/sv_{{svtype}}/batch_{batch}.fa', batch=range(TD_FA_BATCH))
    wildcard_constraints:
        subset='notrsd|sdnotr'
    params:
        cluster_opts='--mem=4096 -t 00:30'
    run:

        fa_pattern = 'temp/td/remap/fa_batch/{subset}/sv_{svtype}/batch_{{batch}}.fa'.format(**wildcards)

        # Read variants
        id_set = set(subset_dtab(
            pd.read_csv(input.tsv_dtab, sep='\t', usecols=('ID', 'IS_TR', 'IS_SD', 'SVTYPE')),
            wildcards.subset, wildcards.svtype
        )['ID'])

        fa_file_name = MERGED_SV_FA.format(vartype='sv', svtype=wildcards.svtype)

        # Read sequences
        record_list = collections.defaultdict(list)

        with gzip.open(fa_file_name, 'rt') as in_file:
            for record in Bio.SeqIO.parse(in_file, 'fasta'):
                if record.id in id_set:
                    record_list[hash(record.id) % TD_FA_BATCH].append(record)

        # Write
        for batch in range(TD_FA_BATCH):
            if len(record_list[batch]) == 0:
                raise RuntimeError(f'FASTA batch {batch} is empty')

            with open(fa_pattern.format(batch=batch), 'wt') as out_file:
                Bio.SeqIO.write(record_list[batch], out_file, 'fasta')

# eshift_td_make_filter_bed
#
# Make BLAST filter region.
rule eshift_td_make_filter_bed:
    output:
        bed='results/td/blast/filter/rmsk_tr_filter.bed.gz'
    params:
        cluster_opts='--mem=2048 -t 00:15'
    shell:
        """zcat {BED_TRF_REGIONS} {BED_RMSK_50} | """
        """grep -Ev '^#' | """
        """cut -f1-3 | """
        """sort -k1,1 -k2,2n | """
        """{bedtools} merge -d 200 | """
        """gzip > {output.bed}"""
