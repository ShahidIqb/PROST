def z_scale(wild,mutant):
    ''' This function gives a z-scale feature by subtracting the z-scale score of wild-residue from mutant residue.
    This function is acquired from iFeature, available at https://github.com/Superzchen/iFeature
    Input: wild-residue, mutant-residue
    Output: Integer vector of five values'''
    amino_dict = {
        'A': [0.24, -2.32, 0.6, -0.14, 1.3],
        'C': [0.84, -1.67, 3.71, 0.18, -2.65],
        'D': [3.98, 0.93, 1.93, -2.46, 0.75],
        'E': [3.11, 0.26, -0.11, -3.04, -0.25],
        'F': [-4.22, 1.94, 1.06, 0.54, -0.62],
        'G': [2.05, 4.06, 0.36, -0.82, -0.38],
        'H': [2.47, 1.95, 0.26, 3.90, 0.09],
        'I': [-3.89, -1.73, -1.71, -0.84, 0.26],
        'K': [2.29, 0.89, -2.49, 1.49, 0.31],
        'L': [-4.28, -1.30, -1.49, -0.72, 0.84],
        'M': [-2.85, -0.22, 0.47, 1.94, -0.98],
        'N': [3.05, 1.60, 1.04, -1.15, 1.61],
        'P': [-1.66, 0.27, 1.84, 0.70, 2.00],
        'Q': [1.75, 0.50, -1.44, -1.34, 0.66],
        'R': [3.52, 2.50, -3.50, 1.99, -0.17],
        'S': [2.39, -1.07, 1.15, -1.39, 0.67],
        'T': [0.75, -2.18, -1.12, -1.46, -0.40],
        'V': [-2.59, -2.64, -1.54, -0.85, -0.02],
        'W': [-4.36, 3.94, 0.59, 3.44, -1.59],
        'Y': [-2.54, 2.44, 0.43, 0.04, -1.47],
    }
    zaamt = np.array(list(amino_dict[mutant]), dtype=np.float)
    zaawt = np.array(list(amino_dict[wild]), dtype=np.float)
    z_scale_desc = zaamt-zaawt
    return list(z_scale_desc)
import numpy as np