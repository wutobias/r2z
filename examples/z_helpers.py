import numpy as np
import pint

ureg = pint.UnitRegistry()

def setLabels(rdmol):
    for atom in rdmol.GetAtoms():
        atom.SetProp('atomLabel',str(atom.GetIdx()))

def is_z_valid(zm, rdmol):
    is_ok        = True
    atm_idx_list = list()
    atm_idx_list.append(zm.ordered_atom_list[0])
    if zm.N_atms != rdmol.GetNumAtoms():
        print(f"Not all atoms in z matrix.")
        is_ok = False
    for z_idx in range(1, zm.N_atms):
        z_row = zm.z[z_idx]
        if None in z_row:
            print(f"Found None at position {z_idx}")
            is_ok = False
        for atm_idx in z_row[1:]:
            if not atm_idx in atm_idx_list:
                print(f"{atm_idx} at position {z_idx} not properly defined")
                is_ok = False
        atm_idx_list.append(z_row[0])
    return is_ok

def calc_rmsd(xyz1, xyz2):
    diff = xyz1.to(ureg.angstrom).magnitude-xyz2.to(ureg.angstrom).magnitude
    rmsd = np.sqrt(np.mean(diff**2))
    return rmsd*ureg.angstrom

def calc_max(xyz1, xyz2):
    diff = xyz1.to(ureg.angstrom).magnitude-xyz2.to(ureg.angstrom).magnitude
    maxdiff = np.max(diff)
    return maxdiff * ureg.angstrom