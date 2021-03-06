"""
Author: Jake Lustig-Yaeger
"""

import os
import sys
import numpy as np
import readsmart as rs
import pandas as pd
from numba import jit

from spectrum_class import Molecule, write_molecular_csv, write_molecular_csv_hr

"""
Create list of molecules and their features
"""
molecules = [
    Molecule("O2", "O_2", [0.2, 0.69, 0.76, 1.27]),
    Molecule("O3", "O_3", [2.51, 0.55, 9.6]),
    Molecule("O2-O2", "O_2-O_2", [0.345, 0.36, 0.38, 0.45, 0.48, 0.53, 0.57, 0.63, 1.065, 1.266]),
    Molecule("CH4", "CH_4", [0.79, 0.89, 1.0, 1.15, 1.380, 1.7, 2.31, 3.3, 7.7]),
    Molecule("CO2", "CO_2", [1.05, 1.21, 1.6, 2.01, 2.730, 4.2, 9.4, 10.4, 15.0]),
    Molecule("CO", "CO", [1.570, 2.35]),
    Molecule("N2-N2", "N_2-N_2", [4.1]),
    Molecule("N2O", "N_2O", [2.11, 2.25, 2.6, 2.67, 2.97, 3.6, 3.9, 4.3, 4.5, 7.9, 17.0]),
    Molecule("H2", "H_2", [0.65, 0.825]),
    Molecule("H2O", "H_2O", [0.651, 0.72, 0.82, 0.94, 1.12, 1.400, 1.85, 2.5, 6.3])
]

def check_molec_uniqueness_for_tableau(molecules = molecules):
    """
    The current Tableau viz fails to correctly label molecular bands that are at
    the same wavelength to within 0.001. This function identifies issues.
    """
    mdic = {}
    # Loop over molecules
    for mol in molecules:
        # Create molecule dictionary
        mdic[mol.name] = mol.bandcenters

    # Check for approximately equal bands between different molecules
    dones = []
    issues = []
    for k1, v1 in mdic.items():
         for k2, v2 in mdic.items():
             tmp = str(sorted("%s%s" %(k1,k2)))
             if k1 != k2 and tmp not in dones:
                 dones.append(tmp)
                 for i, ival in enumerate(v1):
                     for j, jval in enumerate(v2):
                         if np.fabs(ival - jval) < 1e-3:
                             print("WARNING (%s and %s) have (%.3f, %.3f)" %(k1, k2, ival, jval))
                             issues.append((k1, k2, ival, jval))

    if len(issues) == 0:
        print("No issues found.")

    return

if __name__ == "__main__":

    """
    Create test csv file
    """
    check_molec_uniqueness_for_tableau()
    write_molecular_csv(molecules, savename="csv/molecules.csv")
    #write_molecular_csv_hr(molecules, lammin=0.2, lammax=20.0, dwno=1.0, savename="csv/molecules_hr.csv")
