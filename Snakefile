"""
Map insertion sequences back to the reference with BLAST to find patterns of duplications.
"""

import collections
import gzip
import intervaltree
import numpy as np
import os
import pandas as pd
import pysam
import re
import sys

import Bio.SeqIO
import Bio.bgzf

DMAP_DIR = os.path.dirname(os.path.realpath(workflow.snakefile))

sys.path.append(DMAP_DIR)
sys.path.append(os.path.join(DMAP_DIR, 'dep', 'svpop'))
sys.path.append(os.path.join(DMAP_DIR, 'dep', 'svpop', 'dep'))
sys.path.append(os.path.join(DMAP_DIR, 'dep', 'svpop', 'dep', 'ply'))

import libdupmap
import svpoplib


#
# Init
#

# Read and check configuration file

config_filename = libdupmap.pipeline.find_config(config)

if config_filename is None:
    raise RuntimeError('No configuration file found (checked default "config.yaml" and "config.json"')

print(f'Reading config: {config_filename}')

config = None

if config_filename.lower().endswith('.json'):
    import json

    with open(config_filename, 'rt') as config_in:
        config = json.load(config_in)

elif config_filename.lower().endswith('.yaml'):
    import yaml

    with open(config_filename, 'rt') as config_in:
        config = yaml.safe_load(config_in)

else:
    raise RuntimeError(f'Unrecognized configuration file type (expected .json or .yaml): {config_filename}')

if config is None:  # YAML loader returns None if there is nothing in the YAML file
    config = dict()

libdupmap.pipeline.check_config(config)


# Shell prefix
if config['shell_prefix'] != '':
    shell.prefix(config['shell_prefix'])

print(f'shell_prefix: {config["shell_prefix"]}')


#
# Rules
#

localrules: dmap_all

rule dmap_all:
    input:
        bed=expand('results/{sample}/blast/td_sv_ins.bed.gz', sample=libdupmap.pipeline.read_sample_table(config).index)

include: os.path.join(DMAP_DIR, 'rules/input.smk')

include: os.path.join(DMAP_DIR, 'rules/data_ref.smk')
include: os.path.join(DMAP_DIR, 'rules/sample.smk')
include: os.path.join(DMAP_DIR, 'rules/blast.smk')

#include: 'rules/blast.smk'
#include: 'rules/align.smk'
