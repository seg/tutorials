"""
Shim to import data module from the "root" tutorial directory and a function
to return a full path to a file in the image directory.
"""
import sys
import os

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(PARENT_DIR)

# Data module in parent directory
import data

def image_dir(filename):
    """Return the full path to the given filename in the directory where images
    are stored for this tutorial."""
    return os.path.abspath(os.path.join(PARENT_DIR, 'images', filename))
