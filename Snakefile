"""
Map insertion sequences back to the reference with BLAST to find patterns of duplications.
"""

import collections
import gzip
import pandas as pd
import os
import re
import sys

import Bio.SeqIO

sys.path.append(os.path.join(workflow.basedir, 'dep', 'svpop'))
sys.path.append(os.path.join(workflow.basedir, 'dep', 'svpop', 'dep'))
sys.path.append(os.path.join(workflow.basedir, 'dep', 'svpop', 'dep', 'ply'))

import libdupmap
import svpoplib


#
# Init
#

# Read and check configuration file

config_file_name = libdupmap.pipeline.find_config(config)

if config_file_name is not None:
    configfile: config_file_name

libdupmap.pipeline.check_config(config)


# Shell prefix
if config['shell_prefix'] != '':
    shell.prefix(config['shell_prefix'])

print(f'shell_prefix: {config["shell_prefix"]}')


#
# Rules
#

# # eshift_td_blast_filter_basic_del
# #
# # Apply basic filters: DEL
# rule dmap_td_blast_filter_basic_del:
#     input:
#         bed='results/td/blast/full/{subset}/full_td_blast_sv_del.bed.gz',
#         tsv_dtab='results/key/dtab/dtab_sv_del.tsv.gz'
#     output:
#         bed='results/td/blast/filtered/{subset}/filtered_td_blast_sv_del.bed.gz'
#     params:
#         pad=config['td_pad'],
#         pad_min=config['td_pad_min'],
#         chunksize=config['chunksize']
#     run:
#
#         # Read variant BED
#         df_sv = pd.read_csv(config['variant_bed'], sep='\t')
#
#         missing_col = [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN') if col not in df_sv.columns]
#
#         if missing_col:
#             raise RuntimeError(f'Missing required columns in the variant BED file: {", ".join(missing_col)}')
#
#         # Set range where BLAST hits must be found
#         df_sv['PAD'] = df_sv['SVLEN'].apply(lambda val: max([val * params.pad, params.pad_min])).astype(int)
#         df_sv['BLAST_MIN'] = df_sv['POS'] - df_sv['PAD']
#         df_sv['BLAST_MAX'] = df_sv['END'] + df_sv['PAD']
#
#         blast_min = df_sv.set_index('ID')['BLAST_MIN']
#         blast_max = df_sv.set_index('ID')['BLAST_MAX']
#
#         sv_pos = df_sv.set_index('ID')['POS']
#         sv_end = df_sv.set_index('ID')['END']
#
#         # Filter BLAST hits
#         df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=params.chunksize)
#
#         df_list = list()
#
#         for df in df_iter:
#
#             # Filter for BLAST hits between the breakpoints
#             df['BLAST_MIN'] = list(blast_min[df['TIG']])
#             df['BLAST_MAX'] = list(blast_max[df['TIG']])
#
#             df = df.loc[
#                 (df['END'] > df['BLAST_MIN']) &
#                 (df['POS'] < df['BLAST_MAX'])
#             ]
#
#             # Remove direct deletion region hits
#             df['SV_POS'] = list(sv_pos[df['TIG']])
#             df['SV_END'] = list(sv_end[df['TIG']])
#
#             ro = df.apply(lambda row:
#                 svpoplib.variant.reciprocal_overlap(row['POS'], row['END'], row['SV_POS'], row['SV_END']),
#                 axis=1
#             )
#
#             df = df.loc[ro < 0.5]
#
#             # Save
#             if df.shape[0] > 0:
#                 del(df['BLAST_MIN'])
#                 del(df['BLAST_MAX'])
#                 del(df['SV_POS'])
#                 del(df['SV_END'])
#                 df_list.append(df)
#
#         df = pd.concat(df_list, axis=0)
#
#         # Write
#         df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
#
# # eshift_td_blast_filter_basic_ins
# #
# # Apply basic filters: INS
# rule dmap_td_blast_filter_basic_ins:
#     input:
#         bed='blast/blast_full.bed.gz'
#     output:
#         bed='results/td/blast/filtered/{subset}/filtered_td_blast_sv_ins.bed.gz'
#     params:
#         pad=config['td_pad'],
#         pad_min=config['td_pad_min'],
#         chunksize=config['chunksize']
#     run:
#
#         # Read variant BED
#         df_sv = pd.read_csv(config['variant_bed'], sep='\t')
#
#         missing_col = [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN') if col not in df_sv.columns]
#
#         if missing_col:
#             raise RuntimeError(f'Missing required columns in the variant BED file: {", ".join(missing_col)}')
#
#         # Set range where BLAST hits must be found
#         df_sv['PAD'] = df_sv['SVLEN'].apply(lambda val: max([val * params.pad, params.pad_min])).astype(int)
#         df_sv['BLAST_MIN'] = df_sv['POS'] - df_sv['PAD']
#         df_sv['BLAST_MAX'] = df_sv['END'] + df_sv['PAD']
#
#         blast_min = df_sv.set_index('ID')['BLAST_MIN']
#         blast_max = df_sv.set_index('ID')['BLAST_MAX']
#
#         sv_pos = df_sv.set_index('ID')['POS']
#         sv_end = df_sv.set_index('ID')['END']
#
#         # Filter BLAST hits
#         df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=params.chunksize)
#
#         df_list = list()
#
#         for df in df_iter:
#
#             # Filter for BLAST hits between the breakpoints
#             ¡df['BLAST_MIN'] = list(blast_min[df['TIG']])
#             df['BLAST_MAX'] = list(blast_max[df['TIG']])
#
#             df = df.loc[
#                 (df['END'] > df['BLAST_MIN']) &
#                 (df['POS'] < df['BLAST_MAX'])
#             ]
#
#             # Save record
#             if df.shape[0] > 0:
#                 del(df['BLAST_MIN'])
#                 del(df['BLAST_MAX'])
#                 df_list.append(df)
#
#         df = pd.concat(df_list, axis=0)
#
#         # Write
#         df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# Align tasks
#

# Run alignment
rule dmap_align:
    input:
        fa='temp/align/fa/part_{part}.fa',
    output:
        sam=temp('temp/align/sam/align_{part}.sam.gz')
    params:
        aligner=config['aligner'],
        align_params=config['align_params'],
        reference=config['reference']
    threads: 4
    shell:
        """{params.aligner} """
            """{params.align_params} """
            """--secondary=no -a -t {threads} --eqx -Y """
            """{params.reference} {input.fa} | """
            """awk -vOFS="\\t" '($1 !~ /^@/) {{$10 = "*"; $11 = "*"}} {{print}}' | """
            """gzip > {output.sam}"""


#
# BLAST tasks
#

# Merge full BLAST output.
rule dmap_blast_merge_full:
    input:
        bed=expand('temp/blast/part/blast_{part}.bed.gz', part=range(config['partitions']))
    output:
        bed='blast/blast_full.bed.gz'
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
rule dmap_blast:
    input:
        fa='temp/blast/fa/part_{part}.fa',
        bed_filter='data/filter.bed.gz'
    output:
        bed=temp('temp/blast/part/blast_{part}.bed.gz')
    threads: 4
    params:
        temp_dir=lambda wildcards: os.path.join(config['temp_dir'], 'blast_{wildcards.part}'),
        blast_db=config['blast_db'],
        chunksize=config['chunksize'],
        filter_blast=os.path.join(workflow.basedir, 'scripts', 'filter_blast.py')
    shell:
        """blastn -db {params.blast_db} """
            """-num_threads {threads} """
            """-outfmt "6 saccver sstart send qseqid qstart qend sstrand evalue mismatch gaps pident" """
            """-word_size 16 """
            """-perc_identity 95 """
            """-query {input.fa} | """
        """{params.filter_blast} -b {input.bed_filter} -t {params.temp_dir} -c {params.chunksize} -v | """
        """gzip > {output.bed}"""


#
# Initial tasks
#

# Split input sources into partitioned FASTAs
rule dmap_partition_fa:
    input:
        bed='data/variants.bed.gz'
    output:
        fa=expand('temp/{{map_type}}/fa/part_{part}.fa', part=range(config['partitions']))
    wildcard_constraints:
        part=r'\d+',
        map_type='blast|align'
    run:

        out_pattern = f'temp/{wildcards.map_type}/fa/part_{{part}}.fa'

        # Read variants (already tagged by partition)
        df = pd.read_csv(input.bed, sep='\t')

        if wildcards.map_type == 'blast':
            df = df.loc[df['TAG_BLAST']]
        elif wildcards.map_type == 'align':
            df = df.loc[df['TAG_ALIGN']]
        else:
            raise RuntimeError(f'Unknown value for the "map_type" wildcard: {wildcards.map_type}')

        # Write subset FASTA
        for part in range(config['partitions']):

            id_set = set(df.loc[df['PARTITION'] == part, 'ID'])

            try:
                with open(out_pattern.format(part=part), 'wt') as out_file:
                    Bio.SeqIO.write(
                        svpoplib.seq.fa_to_record_iter(config['input'], id_set),
                        out_file,
                        'fasta'
                    )

            except Exception as ex:
                raise RuntimeError(f'Error reading input for partition {part}: {str(ex)}')


# Tag records in variant BED file for ONT or Map.
rule dmap_tag_sv:
    output:
        bed='data/variants.bed.gz'
    run:

        # Read
        df = pd.read_csv(config['variant_bed'], sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'))
        df.set_index('ID', inplace=True, drop=False)
        df.index.name = 'INDEX'

        if len(set(df['ID'])) != df.shape[0]:
            dup_list = sorted([chrom for chrom, count in collections.Counter(df['ID']) if count > 1])
            n_dup = len(dup_list)

            raise RuntimeError('Input variants contain {} duplicate IDs: {}{}'.format(
                len(dup_list), ', '.join(dup_list[:3]), '...' if len(dup_list) > 3 else ''
            ))

        # Check FAI
        df_fai = svpoplib.ref.get_df_fai(config['input_fai'])

        if len(set(df_fai.index)) != df_fai.shape[0]:
            dup_list = sorted([chrom for chrom, count in collections.Counter(df_fai.index) if count > 1])
            n_dup = len(dup_list)

            raise RuntimeError('Input FASTA contains {} duplicate keys: {}{}'.format(
                len(dup_list), ', '.join(dup_list[:3]), '...' if len(dup_list) > 3 else ''
            ))

        missing_set = set(df['ID']) - set(df_fai.index)

        if missing_set:
            raise RuntimeError('Input FASTA is missing {} IDs found in the variant callset: {}{}'.format(
                len(missing_set), ', '.join(sorted(missing_set)[:3]), '...' if len(missing_set) > 3 else ''
            ))

        # Tag variants for BLAST or Map
        blast_range_min = config['blast_range_min']
        blast_range_max = config['blast_range_max']
        align_range_min = config['align_range_min']
        align_range_max = config['align_range_max']

        df['TAG_BLAST'] = True

        df['TAG_BLAST'] = df.apply(lambda row: row['TAG_BLAST'] if row['SVLEN'] >= blast_range_min else False, axis=1)
        if blast_range_max is not None:
            df['TAG_BLAST'] = df.apply(lambda row: row['TAG_BLAST'] if row['SVLEN'] <= blast_range_max else False, axis=1)

        df['TAG_ALIGN'] = True

        df['TAG_ALIGN'] = df.apply(lambda row: row['TAG_ALIGN'] if row['SVLEN'] >= align_range_min else False, axis=1)
        if align_range_max is not None:
            df['TAG_ALIGN'] = df.apply(lambda row: row['TAG_ALIGN'] if row['SVLEN'] <= align_range_max else False, axis=1)

        # Assign partitions
        if config['partitions'] > 1:
            partitions = libdupmap.partition_chrom.partition(df.set_index('ID')['SVLEN'], config['partitions'])
        else:
            partitions = [tuple(sorted(df_fai.index))]

        df['PARTITION'] = pd.Series({
            id: part for part in range(len(partitions)) for id in partitions[part]
        })

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# Make filter BED
rule dmap_make_filter_bed:
    output:
        bed='data/filter.bed.gz'
    threads: 1
    run:

        # Get reference masked locations
        df_mask = libdupmap.ref.masked_fasta_to_bed(
            config['reference'],
            soft=config['ref_soft_mask'],
            dist=config['mask_merge_dist']
        )

        # Merge
        df_mask = libdupmap.ref.merge_regions(
            [df_mask] + config['mask_bed_list'],
            dist=config['mask_merge_dist']
        )

        # Write
        df_mask.to_csv(output.bed, sep='\t', index=False, compression='gzip')
