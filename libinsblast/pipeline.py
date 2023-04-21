"""
Routines for the base pipeline.
"""

import os

import util

def check_config(config):
    """
    Check the configuration object for required elements.
    """

    # Check reference
    if 'reference' not in config.keys():
        raise ValueError('Required key "reference" is not in the configuration')

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

    # Check configuration file name
    config_file_name = config.get('config', None).strip()

    if config_file_name is None or config_file_name.strip() == '':
        return None

    try:
        util.check_regular_file(config_file_name)
    except Exception as ex:
        raise FileNotFoundError(f'Pipeline configuration file not found: {config_file_name}: {str(ex)}')

    if not config_file_name.lower().endswith('.json') or config_file_name.lower().endswith('.json'):
        raise ValueError(f'Expected config file name to end with ".json" or ".yaml": {config_file_name}')

    return config_file_name
