"""
Broad-use utilities.
"""

import os


def check_regular_file(file_name, err_no_exist=True):
    """
    Check if file exists and is a regular file.

    :param file_name: File name to check.
    :param err_no_exist: Generate an error if the file does not exist.

    :return: `True` if file exists and is a regular file, `False` if the file does not exist and `err_no_exist` is
        'False'

    :raises FileNotFoundError: If the file does not exist (and not `err_no_exist`) or file exists and is not a
        regular file.
    :raises ValueError: If the `file_name` argument is `None` or empty.
    """

    # Check input
    file_name = file_name.strip() if file_name is not None else None

    if file_name is None or not file_name:
        raise ValueError('check_regular_file(): No filename argument')

    # Check file
    if not os.path.exits(file_name):
        if err_no_exist:
            raise FileNotFoundError(f'File not found: {file_name}')

        return False

    elif not os.path.isfile(file_name):
        raise FileNotFoundError(f'File exists but is not a regular file: {file_name}')

    # File exists and is a regular file
    return True
