"""
Routines for the base pipeline.
"""

import os

import re

import collections
import libdupmap
import numpy as np
import pandas as pd
import svpoplib


#
# Default vaules
#

DEFAULT_REF_SOFT_MASK = True

DEFAULT_MERGE_DIST = 20

DEFAULT_PARTITIONS = 80

DEFAULT_TD_PAD = 0.1

DEFAULT_TD_PAD_MIN = 100

DEFAULT_CHUNKSIZE = 50000

DEFAULT_BLAST_RANGE = '50:30k'

DEFAULT_ALIGN_RANGE = '50:'

DEFAULT_ALIGNER = 'minimap2'

DEFAULT_ALIGNER_PARAMS = '-x asm5'

DEFAULT_CONFIG_NAME = 'config.yml'

DEFAULT_MIN_SVLEN = 50

DEFAULT_MAX_SVLEN = None

DEFAULT_VCF_FILTER_GT = True

DEFAULT_VCF_PASS = {'PASS', '.'}

DEFAULT_VCF_STRICT_SAMPLE = False

DEFAULT_BLAST_PARAM = '-word_size 16 -perc_identity 95'


#
# Config
#

def check_config(config):
    """
    Check the configuration object for required elements.
    """

    if config is None:
        raise RuntimeError('Configuration object is None')

    # Check reference
    if 'reference' not in config.keys():
        raise ValueError('Required key "reference" is not in the configuration')

    config['reference'] = config['reference'].strip()

    if not config['reference']:
        raise ValueError('Required key "reference" is empty')

    try:
        libdupmap.util.check_regular_file(config['reference'])
    except FileNotFoundError as ex:
        raise FileNotFoundError(f'Required key "reference": {ex}')

    # Check reference FAI
    if 'reference_fai' not in config.keys():
        config['reference_fai'] = config['reference'] + '.fai'

    try:
        libdupmap.util.check_regular_file(config['reference_fai'])
    except FileNotFoundError as ex:
        raise FileNotFoundError(f'Error finding reference FAI: {ex}')


    # Check for sample table
    if 'sample_table' in config.keys():
        if not config['sample_table']:
            raise ValueError('Parameter "sample_table" is empty')

        try:
            libdupmap.util.check_regular_file(config['sample_table'])
        except FileNotFoundError as ex:
            raise FileNotFoundError(f'Parameter "sample_table": {ex}')

    else:
        config['sample_table'] = None

        for filename in ('samples.tsv', 'samples.xlsx'):
            if os.path.exists(filename):
                try:
                    libdupmap.util.check_regular_file(filename)
                except FileNotFoundError as ex:
                    raise FileNotFoundError(f'Default sample table error: {ex}')

                config['sample_table'] = filename

                break

        if config['sample_table'] is None:
            raise RuntimeError('No "sample_table" in config and default sample table filenames were not found (sample.tsv or samples.xlsx)')

    # Set type
    if config['sample_table'].lower().endswith('.tsv') or config['sample_table'].lower().endswith('.tsv.gz'):
        config['sample_table_type'] = 'tsv'
    elif config['sample_table'].lower().endswith('.xlsx') or config['sample_table'].lower().endswith('.xls'):
        config['sample_table_type'] = 'excel'
    else:
        raise RuntimeError(f'Unrecognized file type for parameter "sample_table": Expected ".tsv", ".tsv.gz", ".xlsx", or ".xls": {config["sample_table_type"]}')

    # Reference soft masking
    try:
        config['ref_soft_mask'] = svpoplib.util.as_bool(config.get('ref_soft_mask', DEFAULT_REF_SOFT_MASK))
    except ValueError as ex:
        raise ValueError('Error in config parameter "ref_soft_mask" (expected "true" or "yes")')

    # Mask merge distance
    if 'mask_merge_dist' in config:
        try:
            config['mask_merge_dist'] = libdupmap.util.parse_int(config['mask_merge_dist'])
        except ValueError as ex:
            raise RuntimeError(f'Parameter "mask_merge_dist" is not an integer: {config["mask_merge_dist"]}')

        if config['mask_merge_dist'] < 0:
            raise RuntimeError(f'Parameter "mask_merge_dist" may not be negative: {config["mask_merge_dist"]}')
    else:
        config['mask_merge_dist'] = DEFAULT_MERGE_DIST

    # Mask BED list
    if 'mask_bed_list' not in config:
        config['mask_bed_list'] = list()

    config['mask_bed_list'] = [file_path.strip() for file_path in config['mask_bed_list'] if file_path is not None and file_path.strip != '']

    # Partitions
    if 'partitions' in config:
        try:
            config['partitions'] = libdupmap.util.parse_int(config['partitions'])
        except ValueError as ex:
            raise ValueError(f'Parameter "partitions" must be an integer: {ex}')

    else:
        config['partitions'] = DEFAULT_PARTITIONS

    if config['partitions'] < 1:
        raise RuntimeError(f'Parameter "partitions" must be 1 or greater: {config["partitions"]}')

    # Temporary directory
    if 'temp_dir' not in config:
        config['temp_dir'] = os.path.join(os.getcwd(), 'temp', 'working_temp')
    else:
        config['temp_dir'] = config['temp_dir'].strip()

        if os.path.exists(config['temp_dir']) and not os.path.isdir(config['temp_dir']):
            raise RuntimeError(f'Temporary directory exists and is not a directory: {config["temp_dir"]}')

    if os.path.exists(config['temp_dir']) and os.path.samefile(config['temp_dir'], os.getcwd()):
        raise RuntimeError(f'Temporary directory (config "temp_dir") points to the runtime directory.')

    # BLAST database
    if 'blast_db' not in config:
        raise RuntimeError('Config parameter "blast_db" is missing (PREFIX of the BLAST database, e.g. PREFIX + ".nsq" is the BLAST sequence file)')

    config['blast_db'] = config['blast_db'].strip()

    if not config['blast_db']:
        raise RuntimeError('Config parameter "blast_db" is empty')

    # BLAST parameters
    if 'blast_param' not in config or config['blast_param'].strip() == '':
        config['blast_param'] = DEFAULT_BLAST_PARAM

    # TD pad percent
    if 'td_pad' in config:
        try:
            config['td_pad'] = libdupmap.util.to_prop(config['td_pad'])
        except ValueError as ex:
            raise ValueError(f'Cannot convert "td_pad" parameter to a proportion: {str(ex)}')

        if config['td_pad'] == 0.0:
            raise ValueError(f'Parameter "td_pad_min" must be greater than 0.0')

    else:
        config['td_pad'] = DEFAULT_TD_PAD

    # TD pad min
    if 'td_pad_min' in config:
        try:
            config['td_pad_min'] = libdupmap.util.parse_int(config['td_pad_min'])
        except ValueError as ex:
            raise ValueError(f'Cannot convert "td_pad_min" parameter to an integer: {str(ex)}')

        if config['td_pad_min'] < 0:
            raise ValueError(f'Parameter "td_pad_min" is less than 0: {config["td_pad_min"]}')

    else:
        config['td_pad_min'] = DEFAULT_TD_PAD_MIN

    # Shell prefix
    config['shell_prefix'] = config.get('shell_prefix', '').strip()

    if re.search(r'\bset\ +-euo +pipefail\b', config['shell_prefix']) is None:
        # BASH strict mode
        config['shell_prefix'] = 'set -euo pipefail; ' + config['shell_prefix']

    if not re.search(r';\s*$', config['shell_prefix']):
        # End with semicolon
        config['shell_prefix'] += '; '


    # BLAST SV size range
    if 'blast_range' not in config or not config['blast_range'].split():
        config['blast_range'] = DEFAULT_BLAST_RANGE

    tok = config['blast_range'].split(':')

    if len(tok) != 2:
        raise ValueError(f'Configuration "blast_range" must be two numbers separated by ":": {config["blast_range"]}')

    tok[0] = tok[0].strip()
    tok[1] = tok[1].strip()

    if not tok[0]:
        config['blast_range_min'] = 0
    else:
        try:
            config['blast_range_min'] = libdupmap.util.parse_int(tok[0])
        except ValueError as e:
            raise RuntimeError(f'Error translating min BLAST range to an integer ("blast_range" parameter): {e}')

    if not tok[1]:
        config['blast_range_max'] = None
    else:
        try:
            config['blast_range_max'] = libdupmap.util.parse_int(tok[1])
        except ValueError as e:
            raise RuntimeError(f'Error translating max BLAST range to an integer ("blast_range" parameter): {e}')

    # Align SV size range
    if 'align_range' not in config or not config['align_range'].split():
        config['align_range'] = DEFAULT_ALIGN_RANGE

    tok = config['align_range'].split(':')

    if len(tok) != 2:
        raise ValueError(f'Configuration "align_range" must be two numbers separated by ":": {config["align_range"]}')

    tok[0] = tok[0].strip()
    tok[1] = tok[1].strip()

    if not tok[0]:
        config['align_range_min'] = 0
    else:
        try:
            config['align_range_min'] = libdupmap.util.parse_int(tok[0])
        except ValueError as e:
            raise RuntimeError(f'Error translating min Map range to an integer ("align_range" parameter): {e}')

    if not tok[1]:
        config['align_range_max'] = None
    else:
        try:
            config['align_range_max'] = libdupmap.util.parse_int(tok[1])
        except ValueError as e:
            raise RuntimeError(f'Error translating max Map range to an integer ("align_range" parameter): {e}')

    # Aligner
    if 'aligner' not in config or config['aligner'].strip == '':
        config['aligner'] = DEFAULT_ALIGNER

    if 'align_params' not in config or config['align_params'].strip == '':
        config['align_params'] = DEFAULT_ALIGNER_PARAMS

    # VCF SEQ column
    if 'vcf_seq' in config:
        config['vcf_seq'] = config['vcf_seq'].strip()

        if not config['vcf_seq']:
            config['vcf_seq'] = 'SEQ'
    else:
        config['vcf_seq'] = 'SEQ'

    # VCF filter GT
    try:
        config['vcf_filter_gt'] = svpoplib.util.as_bool(config.get('vcf_filter_gt', DEFAULT_VCF_FILTER_GT))
    except ValueError as ex:
        raise ValueError('Error in config parameter "vcf_filter_gt" (expected "true" or "yes")')

    # VCF PASS
    if 'vcf_pass' in config:
        if issubclass(config['vcf_pass'], str):
            config['vcf_pass'] = {config['vcf_pass']}

            if not str:
                config['vcf_pass'] = set()

        elif issubclass(config['vcf_pass'], list) or issubclass(config['vcf_pass'], set):
            config['vcf_pass'] = set(config['vcf_pass'])

        else:
            raise RuntimeError(f'Unexpected value type for parameter "vcf_pass": Expected string, list, or set: {type(config["vcf_pass"])}')
    else:
        config['vcf_pass'] = DEFAULT_VCF_PASS

    # VCF strict sample
    try:
        config['vcf_strict_sample'] = svpoplib.util.as_bool(config.get('vcf_strict_sample', DEFAULT_VCF_STRICT_SAMPLE))
    except ValueError as ex:
        raise ValueError('Error in config parameter "vcf_strict_sample" (expected "true" or "yes")')


    # Set minimum and maximum SVLEN
    svlen_min = [val for val in [config['blast_range_min'], config['align_range_min']] if val is not None]
    config['svlen_min'] = min(svlen_min) if svlen_min else None

    svlen_max = [val for val in [config['blast_range_max'], config['align_range_max']] if val is not None]
    config['svlen_max'] = max(svlen_max) if len(svlen_max) == 2 else None

    # All OK
    return


def find_config(config=None):
    """
    Find the configuration file to load.

    :param config: Existing config (i.e. from the command-line) used to locate the configuration file
        if key "config" exists. Otherwise, the default `DEFAULT_CONFIG_NAME` is used.

    :raises FileNotFoundError: If the configuration file is not found or is not a regular file.
    :raises ValueError: If the configuration file name does not end in ".json" or ".yaml"
        (not case sensitive).
    """

    # Check arguments
    if config is None:
        config = dict()

    # Get configuration file name from config (if set)
    config_filename = config.get('config', '').strip()

    if config_filename is None or config_filename.strip() == '':
        config_filename = None
        is_default_config = True
    else:
        is_default_config = False

    # Find default config file
    if config_filename is None:
        for filename in ('config.yaml', 'config.json'):
            if os.path.isfile(filename):
                config_filename = filename
                break

    if config_filename is None:
        return None

    # Check file
    try:
        libdupmap.util.check_regular_file(config_filename)
    except Exception as ex:
        raise FileNotFoundError(f'Pipeline configuration file error: {config_filename}: {str(ex)}')

    if not (config_filename.lower().endswith('.json') or config_filename.lower().endswith('.yaml')):
        raise ValueError(f'Expected config file name to end with ".json" or ".yaml": {config_filename}')

    return config_filename


def read_sample_table(config):
    """
    Read and check the sample table. Assumes "sample_table" and "sample_table_type" are set and checked and that
    "sample_table" points to a real file (see `check_config()`)

    :param config: Pipeline configuration file.

    :return: Sample table as a Pandas DataFrame.
    """

    # Read
    if config['sample_table_type'] == 'tsv':
        df_sample = pd.read_csv(config['sample_table'], sep='\t', dtype=str)
    elif config['sample_table_type'] == 'excel':
        df_sample = pd.read_excel(config['sample_table'], dtype=str)
    else:
        raise RuntimeError(f'Program bug: Unknown sample table type: {config["sample_table"]}')

    # Check columns
    missing_cols = [col for col in ('SAMPLE', 'DATA') if col not in df_sample.columns]

    if missing_cols:
        raise RuntimeError(f'Missing column(s) from sample table: {", ".join(missing_cols)}')

    # Check for missing sample names
    is_missing = df_sample['SAMPLE'].apply(lambda val: pd.isnull(val) or val.strip() == '')

    if np.any(is_missing):
        missing_index = list(is_missing[is_missing].index)
        n = len(missing_index)
        missing_str = ', '.join(missing_index[:3]) + ('...' if n > 3 else '')

        raise RuntimeError(f'Found {n} records with missing SAMPLE: {missing_str}')

    # Check for missing data
    is_missing = df_sample['DATA'].apply(lambda val: pd.isnull(val) or val.strip() == '')

    if np.any(is_missing):
        missing_index = [val + 1 for val in is_missing[is_missing].index]
        n = len(missing_index)
        missing_str = ', '.join(missing_index[:3]) + ('...' if n > 3 else '')

        raise RuntimeError(f'Found {n} records with missing DATA: {missing_str}')

    # Check for duplicate sample names
    if len(set(df_sample['SAMPLE'])) != df_sample.shape[0]:
        dup_list = sorted([sample for sample, count in collections.Counter(df_sample['SAMPLE']).items() if count > 1])
        n = len(dup_list)
        dup_str = ', '.join(dup_list[:3]) + ('...' if n > 3 else '')

        raise RuntimeError(f'Found {n} duplicate samle names: {dup_str}')

    # Return table
    return df_sample.set_index('SAMPLE')


def get_sample_entry(sample_name, config):
    """
    Get an record from the sample config as a dictionary.

    :param sample_name: Sample name.
    :param config: Pipeline config.

    :return: Sample entry dictionary.
    """

    # Get sample row
    df_sample = read_sample_table(config)

    if sample_name not in df_sample.index:
        raise RuntimeError(f'Sample is not in the sample table: {sample_name}')

    sample_row = df_sample.loc[sample_name]

    sample_data_lower = sample_row['DATA'].lower()

    # BED/FASTA
    if ';' in sample_row['DATA'] or sample_data_lower.endswith('.bed.gz') or sample_data_lower.endswith('.bed'):
        bed_filename = None
        fa_filename = None

        for filename in sample_row['DATA'].split(';'):
            filename = filename.strip()

            if not filename:
                continue

            filename_lower = filename.lower()

            if filename_lower.endswith('.bed') or filename_lower.endswith('.bed.gz'):
                if bed_filename is None:
                    bed_filename = filename
                else:
                    raise RuntimeError(f'Multiple BED files for sample {sample_name}: {bed_filename}, {filename}')

            elif filename_lower.endswith('.fa') or filename_lower.endswith('.fa.gz') or \
                filename_lower.endswith('.fasta') or filename_lower.endswith('.fasta.gz') or \
                filename_lower.endswith('.fna') or filename_lower.endswith('.fna.gz'):

                if fa_filename is None:
                    fa_filename = filename
                else:
                    raise RuntimeError(f'Multiple FASTA files for sample {sample_name}: {fa_filename}, {filename}')

            else:
                raise RuntimeError(f'Unexpected file extension for sample {sample_name}: Expected BED or FASTA files: {filename}')

        # Check BED filename
        if bed_filename is None:
            raise RuntimeError(f'Missing BED filename for sample {sample_name}')

        try:
            libdupmap.util.check_regular_file(bed_filename)
        except FileNotFoundError as ex:
            raise FileNotFoundError(f'BED file error for sample {sample_name}: {ex}')

        # Get BED columns
        with svpoplib.seq.PlainOrGzReader(bed_filename) as in_file:
            line = next(in_file).strip()

        if line:
            bed_cols = line.split('\t')
        else:
            raise RuntimeError(f'Could not read column headers from {bed_filename}')

        bed_has_seq = 'SEQ' in bed_cols

        # Check FASTA filename
        if fa_filename is not None:
            try:
                libdupmap.util.check_regular_file(fa_filename)
            except FileNotFoundError as ex:
                raise FileNotFoundError(f'FASTA file error for sample {sample_name}: {ex}')

            seq_source = 'fasta'

        elif not bed_has_seq:

            # Try to find filename based on SV-Pop directory structures
            fa_filename = os.path.join(
                os.path.dirname(bed_filename),
                'fa',
                re.sub('\.bed(\.gz)?$', '', os.path.basename(bed_filename)) + '.fa.gz'
            )

            try:
                libdupmap.util.check_regular_file(fa_filename)
            except FileNotFoundError as ex:
                raise FileNotFoundError(f'Error finding FASTA (no SEQ column in BED) file assuming SV-Pop relative path from BED file for sample {sample_name}: {ex}')

            seq_source = 'fasta'

        else:
            seq_source = 'bed'

        # Set entry
        return {
            'sample': sample_name,
            'data': sample_row['DATA'],
            'type': 'bed',
            'bed_file': bed_filename,
            'fa_file': fa_filename,
            'bed_cols': bed_cols,
            'seq_source': seq_source
        }

    # VCF
    if sample_data_lower.endswith('.vcf.gz') or sample_data_lower.endswith('.vcf') or sample_data_lower.endswith('.bcf'):

        try:
            libdupmap.util.check_regular_file(sample_row['DATA'])
        except FileNotFoundError as ex:
            raise FileNotFoundError(f'FASTA file error for sample {sample_name}: {ex}')

        return {
            'sample': sample_name,
            'data': sample_row['DATA'],
            'type': 'vcf',
            'vcf_file': sample_row['DATA']
        }

    # Unrecognized type
    raise RuntimeError(f'Unrecognized type in DATA for sample {sample_name}: Expected a BED or VCF file: {sample_row["DATA"]}')
