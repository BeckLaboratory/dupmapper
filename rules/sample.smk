# Sample data preparation

# Split input sources into partitioned FASTAs
rule dmap_sample_partition_fa:
    input:
        bed='temp/{sample}/input/sv_ins.bed.gz',
        fa='temp/{sample}/input/sv_ins.fa.gz'
    output:
        fa=expand('temp/{{sample}}/input/fa_part/part_{part}.fa', part=range(config['partitions']))
    wildcard_constraints:
        part=r'\d+'
    threads: 1
    run:

        out_pattern = f'temp/{wildcards.sample}/input/fa_part/part_{{part}}.fa'

        # Read variants (already tagged by partition)
        df = pd.read_csv(input.bed, sep='\t')

        # Write subset FASTA
        for part in range(config['partitions']):

            id_set = set(df.loc[df['PARTITION'] == part, 'ID'])

            if id_set:
                try:
                    with open(out_pattern.format(part=part), 'wt') as out_file:
                        Bio.SeqIO.write(
                            svpoplib.seq.fa_to_record_iter(input.fa, id_set),
                            out_file,
                            'fasta'
                        )

                except Exception as ex:
                    raise RuntimeError(f'Error reading input for partition {part}: {str(ex)}', ex)

            else:
                # Empty file if no variants
                with open(out_pattern.format(part=part), 'wt') as out_file:
                    pass


# Tag records in variant BED file for BLAST or Map.
rule dmap_sample_tag_sv:
    input:
        bed=dmap_sample_input_bed,
        fa=dmap_sample_input_fa
    output:
        bed='temp/{sample}/input/sv_ins.bed.gz',
        fa='temp/{sample}/input/sv_ins.fa.gz'
    threads: 1
    run:

        sample_entry = libdupmap.pipeline.get_sample_entry(wildcards.sample, config)

        # Read
        usecols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

        if not input.fa:
            usecols += ['SEQ']

        df = pd.read_csv(input.bed, sep='\t', usecols=usecols, dtype={'SEQ': str, '#CHROM': str})
        df.set_index('ID', inplace=True, drop=False)
        df.index.name = 'INDEX'

        # Subset for insertions
        df = df.loc[df['SVTYPE'] == 'INS'].copy()

        # Check for duplicate IDs
        if len(set(df['ID'])) != df.shape[0]:
            dup_list = sorted([chrom for chrom, count in collections.Counter(df['ID']).items() if count > 1])
            n_dup = len(dup_list)

            raise RuntimeError('Input variants contain {} duplicate IDs: {}{}'.format(
                len(dup_list), ', '.join(dup_list[:3]), '...' if len(dup_list) > 3 else ''
            ))

        # Tag variants for BLAST or Map
        blast_range_min = config['blast_range_min']
        blast_range_max = config['blast_range_max']
        align_range_min = config['align_range_min']
        align_range_max = config['align_range_max']

        if blast_range_min is None:
            blast_range_min = 0

        if blast_range_max is None:
            blast_range_max = np.max(df['SVLEN'])

        if align_range_min is None:
            align_range_min = 0

        if align_range_max is None:
            align_range_max = np.max(df['SVLEN'])

        df['TAG_BLAST'] = (df['SVLEN'] >= blast_range_min) & (df['SVLEN'] <= blast_range_max)
        df['TAG_ALIGN'] = (df['SVLEN'] >= align_range_min) & (df['SVLEN'] <= align_range_max)

        df = df.loc[df['TAG_BLAST'] | df['TAG_ALIGN']].copy()

        # Write FASTA
        if input.fa:
            with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                Bio.SeqIO.write(
                    svpoplib.seq.fa_to_record_iter(input.fa, record_set=set(df['ID']), require_all=True, input_format='fasta'),
                    fa_file, 'fasta'
                )

        else:
            # Write sequences

            if 'SEQ' not in df.columns:
                raise RuntimeError('No variant sequence FASTA or SEQ in input BED file')

            if np.any(pd.isnull(df['SEQ'])):
                n_null = np.sum(pd.isnull(df['SEQ']))

                raise RuntimeError(f'Found {n_null} records in BED with a null SEQ')

            with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                Bio.SeqIO.write(
                    svpoplib.seq.bed_to_seqrecord_iter(df),
                    fa_file, 'fasta'
                )

            del(df['SEQ'])

        # Assign partitions
        if config['partitions'] > 1:
            partitions = libdupmap.partition_chrom.partition(df.set_index('ID')['SVLEN'], config['partitions'])
        else:
            partitions = [tuple(sorted(df['ID']))]

        df['PARTITION'] = pd.Series({
            id: part for part in range(len(partitions)) for id in partitions[part]
        })

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# VCF to BED
#

rule dmap_sample_vcf_to_bed:
    input:
        vcf=dmap_sample_input_vcf
    output:
        bed='temp/{sample}/vcf_input/sv_ins.bed.gz',
        fa='temp/{sample}/vcf_input/sv_ins.fa.gz',
        fai='temp/{sample}/vcf_input/sv_ins.fa.gz.fai'
    params:
        chunk_size=20000  # DataFrame chunk size
    threads: 12
    run:

        seq_col_name = config['vcf_seq']

        # Build query string
        has_info_seq = False

        with svpoplib.seq.PlainOrGzReader(input.vcf) as in_file:
            for line in in_file:
                line = line.strip()

                if not line:
                    continue

                if not line.startswith('#'):
                    break

                if line.startswith(f'##INFO=<ID={seq_col_name},'):
                    has_info_seq
                    break

        query_string = r'"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER'

        if has_info_seq:
            query_string += r'\tINFO/' + seq_col_name

        has_gt = config['vcf_filter_gt'] or config['vcf_strict_sample']

        if has_gt:
            query_string += r'[\t%GT]'


        query_string += r'\n"'

        # Run bcftools
        table_file_name = os.path.join(os.path.dirname(output.bed), 'bcftools_out.tsv.gz')

        shell(
            """bcftools query -H -f{query_string} {input.vcf} | """
            """gzip """
            """> {table_file_name}"""
        )

        # Parse
        chrom_set = set(svpoplib.ref.get_df_fai(config['reference_fai']).index)

        write_header = True

        retain_cols = [
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'QUAL', 'FILTER'
        ] + (
            ['GT'] if has_gt else []
        ) + [
            'VCF_POS', 'VCF_IDX'
        ]

        try:
            with gzip.open(output.bed, 'wt') as bed_file:
                with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                    for df in svpoplib.variant.vcf_tsv_to_bed(
                        tsv_in=table_file_name,
                        sample=wildcards.sample,
                        bed_file=None,
                        chunk_size=params.chunk_size,
                        threads=threads,
                        filter_pass_set=config['vcf_pass'],
                        filter_gt=config['vcf_filter_gt'],
                        strict_sample=config['vcf_strict_sample']
                    ):

                        # Drop chromosomes not in this reference (e.g. Called against and ALT reference with SV-Pop using a no-ALT reference)
                        df = df.loc[df['#CHROM'].isin(chrom_set)]

                        df = df.loc[(df['SVTYPE'] == 'INS')]

                        if config['svlen_min'] is not None:
                            df = df.loc[df['SVLEN'] >= config['svlen_min']]

                        if config['svlen_max'] is not None:
                            df = df.loc[df['SVLEN'] <= config['svlen_max']]

                        if df.shape[0] > 0:

                            # Write FASTA
                            Bio.SeqIO.write(
                                svpoplib.seq.bed_to_seqrecord_iter(df),
                                fa_file, 'fasta'
                            )

                            # Write BED
                            df[retain_cols].to_csv(bed_file, sep='\t', index=False, header=write_header)
                            bed_file.flush()

                            write_header = False

                    # Write empty BED if no records were found
                    if write_header:
                        raise RuntimeError(f'No records found in VCF: {input.vcf}')
        finally:
            if os.path.isfile(table_file_name):
                os.remove(table_file_name)

        pysam.faidx(output.fa)
