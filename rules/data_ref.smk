# Reference preparation

# Make filter BED
rule dmap_ref_filter_bed:
    output:
        bed='data/ref/filter.bed.gz'
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
