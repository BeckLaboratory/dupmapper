"""
Broad-use utilities.
"""

import os
import numpy as np


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
    if not os.path.exists(file_name):
        if err_no_exist:
            raise FileNotFoundError(f'File not found: {file_name}')

        return False

    elif not os.path.isfile(file_name):
        raise FileNotFoundError(f'File exists but is not a regular file: {file_name}')

    # File exists and is a regular file
    return True


def to_prop(val):
    """
    Convert a value to a floating-point proportion (between 0.0 and 1.0 inclusive).

    :param val: Value to convert to a proportion. If float and in bounds, return the same value. Converted to a
        floating point value if it is a string. If the string ends with "%", then the value is converted to a float
        and divided by 100 to convert a percentage to a proportion.

    :return: Proportion as a Python float.
    """

    if isinstance(val, float) or isinstance(val, np.float):
        float_val = val

    if isinstance(val, np.floating):
        float_val = int(val)

    elif isinstance(val, str):
        val = val.strip()

        if not val:
            raise ValueError('Proportion string is empty')

        if val.endswith('%'):
            try:
                float_val = float(val[:-1]) / 100
            except ValueError as ex:
                raise ValueError(f'Cannot convert percentage (expected "X%" or "X.Y%"): {str(ex)}')

        else:
            try:
                float_val = float(val)
            except ValueError as ex:
                raise ValueError(f'Cannot convert proportion (expected floating point): {str(ex)}')

    else:
        raise ValueError('Proportion is not floating point or a percentage string')

    # Check value
    if not 0.0 <= float_val <= 1.0:
        raise ValueError(f'Proportion is out of bounds (expected [0.0, 1.0]: {val}')

    return float_val


def parse_int(str_val):
    """
    Get an integer value from a string with optional "k", "m", and "g" suffixes (not case sensitive) or scientific
    notation (e.g. 1.2e3).

    :param str_val: String value.

    :return: Integer value.

    :throws ValueError: If the string does not appear to be an integer.
    """

    str_val_lower = str_val.lower()

    if str_val_lower.endswith('k'):
        multiplier = int(1e3)
        str_val = str_val[:-1]

    elif str_val_lower.endswith('m'):
        multiplier = int(1e6)
        str_val = str_val[:-1]

    elif str_val_lower.endswith('g'):
        multiplier = int(1e9)
        str_val = str_val[:-1]

    else:
        multiplier = 1

    return int(float(str_val) * multiplier)
