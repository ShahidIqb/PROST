from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import numpy as np
from itertools import combinations
from collections import Counter

def dssp_features(file_in, chain, position):
    """This function extracts RSA, phi, and psi angle for mutation position."""
    p = PDBParser(QUIET=True)
    structure = p.get_structure('xxx', file_in)
    model = structure[0]
    dssp = DSSP(model, file_in, dssp='mkdssp')
    #print(dssp[(chain, (' ', position, ' '))][3:6])
    return list(dssp[(chain, (' ', position, ' '))][3:6])


def translate_3aa1(three_letter):
    if len(three_letter) > 3:
        three_letter = three_letter[1:4]
    trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
             'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
             'TYR': 'Y', 'VAL': 'V'}
    return (trans[three_letter])


def distance_euc(a, b):
    return np.linalg.norm(a-b)

def read_pdbBlock(path_data, pdb_id, chain, residue, residue_position, min_distance, max_distance):
    """Read a block of protein in the given chain. The geometric mean of the wild residue at given position is calculate and all
    atoms are considered from minimum distance to maximum distance.
    Input: pdb_id, residue, residue_position, min_distance, max_distance
    Output: A block of PDB"""
    #print(path_data, pdb_id, chain, residue, residue_position, min_distance, max_distance)
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0
    d_count = 0
    for line in open(path_data+pdb_id+'.pdb'):
        list_elem = line.split()
        if list_elem[0]=='ATOM' and list_elem[5]==residue_position and translate_3aa1(list_elem[3])==residue and list_elem[4]==chain:
            x_cord =''
            y_cord=''
            z_cord=''
            if len(list_elem[6]) > 8:
                xy_comb = list_elem[6].split('-')
                if xy_comb[0]=='':
                    x_cord = '-'+xy_comb[1]
                    y_cord = '-' + xy_comb[2]
                else:
                    x_cord = xy_comb[0]
                    y_cord = '-'+xy_comb[1]
            else:
                x_cord = list_elem[6]
            if len(list_elem[7])>8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0]=='':
                    y_cord = '-'+yz_comb[1]
                    z_cord = '-' + yz_comb[2]
                else:
                    y_cord = yz_comb[0]
                    z_cord = '-'+yz_comb[1]
            else:
                y_cord = list_elem[7]
            if z_cord=='':
                z_cord = list_elem[8]
            x_sum += float(x_cord)
            y_sum += float(y_cord)
            z_sum += float(z_cord)
            #print(list_elem[3], list_elem[4], list_elem[5], x_cord, y_cord, z_cord)
            d_count += 1
        if d_count > 0 and (line.split()[0]=='TER' or line.split()[0]=='END'):
            break
    xyz_mean = np.array((x_sum/d_count, y_sum/d_count, z_sum/d_count))

    #print(x_sum, y_sum, z_sum)

    pdb_block, pdb_block_list = extract_residueaswhole_within_tresh(path_data,pdb_id, chain, xyz_mean, min_distance, max_distance)
    return pdb_block_list

def extract_residueaswhole_within_tresh(path_data, pdb_id, chain, xyz_mean, min_distance, max_distance):
    """This fucntion reads all the atoms of the residue when CA is in the threshold"""
    residue_details = set()
    for line in open(path_data + pdb_id + '.pdb'):
        line_elem = line.split()
        list_elem = line_elem
        if (line_elem[0] == 'ATOM' and line_elem[4]==chain and line_elem[2]=='CB') or (line_elem[0] == 'ATOM' and line_elem[4]==chain and line_elem[3]=='GLY' and line_elem[2]=='CA'):
            x_cord = 0.0
            y_cord = 0.0
            z_cord = 0.0
            if len(list_elem[6]) > 8:
                xy_comb = list_elem[6].split('-')
                if xy_comb[0] == '':
                    x_cord = float('-' + xy_comb[1])
                    y_cord = float('-' + xy_comb[2])
                else:
                    x_cord = float(xy_comb[0])
                    y_cord = float('-' + xy_comb[1])
            else:
                x_cord = float(line_elem[6])
            if len(list_elem[7]) > 8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0] == '':
                    y_cord = float('-' + yz_comb[1])
                    z_cord = float('-' + yz_comb[2])
                else:
                    y_cord = float(yz_comb[0])
                    z_cord = float('-' + yz_comb[1])
            else:
                y_cord = float(line_elem[7])
            if z_cord==0.0:
                z_cord = float(line_elem[8])
            #print(x_cord, y_cord, z_cord)
            xyz = np.array((x_cord, y_cord, z_cord))
            if min_distance < distance_euc(xyz, xyz_mean) <= max_distance:
                residue_details.add(line_elem[4]+line_elem[5])
        if len(residue_details) > 0 and (line.split()[0]=='TER' or line.split()[0]=='END'):
            break
    #print(residue_details)
    pdb_block_list = []
    pdb_block = ''
    for linee in open(path_data + pdb_id + '.pdb'):
        line_elem = linee.split()
        if line_elem[0] =='ATOM' and line_elem[4]==chain:
            if line_elem[4]+line_elem[5] in residue_details:
                pdb_block = pdb_block+linee
                pdb_block_list.append(linee)
        if len(pdb_block) > 0 and (linee.split()[0]=='TER' or linee.split()[0]=='END'):
            break
    return pdb_block, pdb_block_list

def atom_coordinates(PDB_block):
    '''Coordinates for polar and hydrophobic Ca atoms from the PDB_Block'''
    #print(PDB_block)
    #residue groups from https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
    polar_group = ('R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T', 'Y')
    nonpolar_group = ('A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V')
    hydrophobic_group = ('A', 'C', 'I', 'L', 'M', 'F', 'W', 'V')
    label_coordinates = []
    for line in PDB_block:
        list_elem = line.split()
        if list_elem[2] == 'CA':
            x_cord = 0.0
            y_cord = 0.0
            z_cord = 0.0
            if len(list_elem[6]) > 8:
                xy_comb = list_elem[6].split('-')
                if xy_comb[0] == '':
                    x_cord = float('-' + xy_comb[1])
                    y_cord = float('-' + xy_comb[2])
                else:
                    x_cord = float(xy_comb[0])
                    y_cord = float('-' + xy_comb[1])
            else:
                x_cord = float(list_elem[6])
            if len(list_elem[7]) > 8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0] == '':
                    y_cord = float('-' + yz_comb[1])
                    z_cord = float('-' + yz_comb[2])
                else:
                    y_cord = float(yz_comb[0])
                    z_cord = float('-' + yz_comb[1])
            else:
                y_cord = float(list_elem[7])
            if z_cord == 0.0:
                z_cord = float(list_elem[8])
            if translate_3aa1(list_elem[3]) in hydrophobic_group:
                label_coordinates.append(('H', (float(x_cord), float(y_cord), float(z_cord))))
            if translate_3aa1(list_elem[3]) in polar_group:
                label_coordinates.append(('D', (float(x_cord), float(y_cord), float(z_cord))))
            if translate_3aa1(list_elem[3]) in nonpolar_group:
                label_coordinates.append(('A', (float(x_cord), float(y_cord), float(z_cord))))
    return label_coordinates

def distance_xyz(p1, p2):
    """Euclidean distance between labels with two points."""
    label1, x1y1z1 = p1
    x1, y1, z1 = x1y1z1
    label2, x2y2z2 = p2
    x2, y2, z2 = x2y2z2
    return (label1+'_'+label2, distance_euc(np.array([x1, y1, z1]), np.array([x2, y2, z2]))) #hypot(x2 - x1, y2 - y1, z2 - z1))

def getFrequency(distMatrix, d_step):
    """Frquency of atom classes with cutoff value"""
    atom_classes3 = {'A_A':0, 'A_D':0, 'D_A':0, 'A_H':1, 'H_A':1,
                    'D_D':0, 'D_H':1, 'H_D':1,
                    'H_H':2}
    cmbnation = 3
    scan_matrix = [0] * cmbnation * int(8/d_step)
    i = 0
    for dist in range(2, 8+1, d_step):
        for rec in distMatrix:
            if rec[1] <= dist and rec[0] in atom_classes3:
                scan_matrix[atom_classes3.get(rec[0])+i*cmbnation] += 1
        i+=1
    return scan_matrix


def pharmacophore_count(atom_xyz):
    """This function returns the pharmacophore count"""
    p_count = {}
    for atom in atom_xyz:
        if atom[0] in p_count:
            p_count[atom[0]] += 1
        else:
            p_count[atom[0]] = 1
    return p_count


def pahrmaco_sign(data_path, pdb_id, chain, wild_res, mut_pos):
    """Generates the mCSM signature in the vecinity of mutation site."""
    pdb_block_list = read_pdbBlock(data_path, pdb_id, chain, wild_res, str(mut_pos), 0, 8)
    pharm_coords_simple = atom_coordinates(PDB_block=pdb_block_list)
    distances_simple = []
    for combo in combinations(pharm_coords_simple, 2):
        distances_simple.append(distance_xyz(*combo))
    mCSMFreq_simple = getFrequency(distances_simple, 2)
    pharmaco_count = pharmacophore_count(pharm_coords_simple)
    return mCSMFreq_simple, pharmaco_count

