"""
Map insertion sequences back to the reference with BLAST to find patterns of duplications.
"""

import collections
import gzip
import intervaltree
import numpy as np
import os
import pandas as pd
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


##
## Rules
##

include: 'rules/init.sm'
include: 'rules/blast.sm'
include: 'rules/align.sm'
