# Rule input functions

def dmap_sample_input_bed(wildcards):
    """
    Get the path to the input BED file (direct input or parsed from VCF input).

    :param wildcards: Rule wildcards.

    :return: Path to input BED file.
    """
    sample_entry = libdupmap.pipeline.get_sample_entry(wildcards.sample, config)

    if sample_entry['type'] == 'bed':
        return sample_entry['bed_file'].format(**wildcards)

    elif sample_entry['type'] == 'vcf':
        return 'temp/{sample}/vcf_input/sv_ins.bed.gz'.format(**wildcards)

    raise RuntimeError(f'Program bug: Error getting BED input: Unrecognized sample entry type {sample_entry["type"]}')


def dmap_sample_input_fa(wildcards):
    """
    Get the path to the input FASTA file (direct input or parsed from VCF input).

    :param wildcards: Rule wildcards.

    :return: Path to input FASTA file.
    """
    sample_entry = libdupmap.pipeline.get_sample_entry(wildcards.sample, config)

    if sample_entry['type'] == 'bed':
        if sample_entry['seq_source'] == 'fasta':
            return sample_entry['fa_file'].format(**wildcards)
        else:
            return []

    elif sample_entry['type'] == 'vcf':
        return 'temp/{sample}/vcf_input/sv_ins.fa.gz'.format(**wildcards)

    raise RuntimeError(f'Program bug: Error getting FASTA input: Unrecognized sample entry type {sample_entry["type"]}')


def dmap_sample_input_vcf(wildcards):
    """
    Get the path to the input VCF file.

    :param wildcards: Rule wildcards.

    :return: Path to input VCF file.
    """

    sample_entry = libdupmap.pipeline.get_sample_entry(wildcards.sample, config)

    if 'vcf_file' not in sample_entry:
        raise RuntimeError('No VCF input for sample')

    return sample_entry['vcf_file']