import os
import numpy as np

basedir = os.path.dirname(__file__)
seis_data_name = os.path.join(basedir, 'penobscot_subset.bin')
seis_header_name = os.path.join(basedir, 'penobscot_subset.hdr')
horizon_name = os.path.join(basedir, 'Hor_B.txt')

def load_seismic():
    """Read seismic data into memory as a 3D numpy array."""
    # Seismic stored as a "raw" binary array instead of .npy for compatibility
    # with other languages.
    with open(seis_header_name, 'r') as infile:
        next(infile)
        shape = [int(next(infile).split()[1]) for _ in range(3)]
    data = np.fromfile(seis_data_name, dtype=np.int16)
    return data.reshape(shape).astype(float)

def load_horizon():
    """Read horizon as 2D array of z-indices into seismic cube."""
    return np.loadtxt(horizon_name)

__all__ = ['load_seismic', 'load_horizon']
