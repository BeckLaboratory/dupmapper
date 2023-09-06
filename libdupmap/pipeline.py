"""
Routines for the base pipeline.
"""

import os

import re

import libdupmap
import svpoplib


#
# Default vaules
#

DEFAULT_REF_SOFT_MASK = True

DEFAULT_MERGE_DIST = 20

DEFAULT_PARTITIONS = 80

DEFAULT_BLAST_THREADS = 12

DEFAULT_TD_PAD = 0.1

DEFAULT_TD_PAD_MIN = 100

DEFAULT_CHUNKSIZE = 50000

DEFAULT_BLAST_RANGE = ':30k'

DEFAULT_ALIGN_RANGE = ':'

DEFAULT_ALIGNER = 'minimap2'

DEFAULT_ALIGNER_PARAMS = '-x asm5'


#
# Config
#

def check_config(config):
    """
    Check the configuration object for required elements.
    """

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

    # config['reference_fai'] = config['reference'] + '.fai'
    #
    # try:
    #     libdupmap.util.check_regular_file(config['reference_fai'])
    # except FileNotFoundError as ex:
    #     raise FileNotFoundError(f'Required key "reference" FAI index: {ex}')

    # Input sequence files
    if 'input' not in config:
        raise ValueError('Required key "input" is not in the configuration')

    config['input'] = config['input'].strip()

    if not config['input']:
        raise RuntimeError('No files found in "input" configuration option')

    if re.search(re.compile(r'.*\.fa|fasta(\.gz)?$', re.I), config['input']) is None:
        raise FileNotFoundError(f'Input file ("input" parameter) does not look like a FASTA file name (.fa or .fasta, optional .gz): {config["input"]}')

    try:
        libdupmap.util.check_regular_file(config['input'])
    except FileNotFoundError as ex:
        raise FileNotFoundError(f'Input file ("input" parameter) does not exist or is not a regular file: {config["input"]}')

    config['input_fai'] = config['input'] + '.fai'

    try:
        libdupmap.util.check_regular_file(config['input_fai'])
    except FileNotFoundError as ex:
        raise FileNotFoundError(f'Missing ".fai" index for input file ("input" parameter): {config["input"]}')

    if config['input'].endswith('.gz'):
        try:
            libdupmap.util.check_regular_file(config['input'] + '.gzi')
        except FileNotFoundError as ex:
            raise FileNotFoundError(f'Missing ".gzi" index for gzipped input file ("input" parameter): {config["input"]}')

    # Reference soft masking
    try:
        config['ref_soft_mask'] = svpoplib.util.as_bool(config.get('ref_soft_mask', DEFAULT_REF_SOFT_MASK))
    except ValueError as ex:
        raise ValueError('Error in config parameter "ref_soft_mask" (expected "true" or "yes")')

    # Mask merge distance
    if 'mask_merge_dist' in config:
        try:
            config['mask_merge_dist'] = int(config['mask_merge_dist'])
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
            config['partitions'] = int(config['partitions'])
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

    # BLAST threads
    if 'blast_threads' in config:
        try:
            config['blast_threads'] = int(config['blast_threads'])
        except ValueError as ex:
            raise RuntimeError(f'Config parameter "blast_threads" is not an integer: {str(ex)}')

        if config['blast_threads'] < 1:
            raise RuntimeError(f'Config parameter "blast_threads" must be at least 1: {config["blast_threads"]}')
    else:
        config['blast_threads'] = DEFAULT_BLAST_THREADS

    # BLAST database
    if 'blast_db' not in config:
        raise RuntimeError('Config parameter "blast_db" is missing (PREFIX of the BLAST database, e.g. PREFIX + ".nsq" is the BLAST sequence file)')

    config['blast_db'] = config['blast_db'].strip()

    if not config['blast_db']:
        raise RuntimeError('Config parameter "blast_db" is empty')

    # Variant BED file
    if 'variant_bed' not in config:
        raise RuntimeError('Config parameter "variant_bed" is missing (variants in BED format with columns: #CHROM, POS, END, SVTYPE, SVLEN)')

    config['variant_bed'] = config['variant_bed'].strip()

    if not config['variant_bed']:
        raise RuntimeError('Config parameter "variant_bed" is empty')

    try:
        libdupmap.util.check_regular_file(config['variant_bed'])
    except FileNotFoundError as ex:
        raise FileNotFoundError(f'Missing variant BED ("variant_bed" parameter): {str(ex)}')

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
            config['td_pad_min'] = int(config['td_pad_min'])
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

    # All OK
    return


# def expand_input(input_list, delim=';', strip=True):
#     """
#     Read a list of values and generate a flat list of input files. The input list may be a single path (str),
#         multiple paths separated by `delim`, or a collection of paths (list, dict, tuple, etc). Each element that
#         is not type `str` is treated as a collection. Each collection may contain collections.
#
#     :param input_list: A collection of input files (see description above).
#     :param delim: Split strings on this delimiter if not `None` where each element in the split string is a separate
#         input path.
#     :param strip: Strip whitespace from filenames.
#
#     :return: A flat list of input files separated from the collections and split on
#     """
#
#     exp_list = list()
#
#     if isinstance(input_list, str):
#         input_list = [input_list]
#     else:
#         input_list = list(input_list).copy()
#
#     while len(input_list) > 0:
#         input_item = input_list.pop(0)
#
#         if isinstance(input_item, str):
#             if delim is not None:
#                 for subitem in input_item.split(delim):
#                     if strip:
#                         subitem = subitem.strip()
#
#                     if subitem:
#                         exp_list.append(subitem)
#                 else:
#                     exp_list.append(input_item)
#
#         else:
#             for val in input_item[::-1]:
#                 input_list.insert(0, val)
#
#     return exp_list


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

    # Check configuration file name
    config_file_name = config.get('config', '').strip()

    if config_file_name is None or config_file_name.strip() == '':
        return None

    try:
        libdupmap.util.check_regular_file(config_file_name)
    except Exception as ex:
        raise FileNotFoundError(f'Pipeline configuration file not found: {config_file_name}: {str(ex)}')

    if not config_file_name.lower().endswith('.json') or config_file_name.lower().endswith('.yaml'):
        raise ValueError(f'Expected config file name to end with ".json" or ".yaml": {config_file_name}')

    return config_file_name
