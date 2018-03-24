#!/usr/bin/python3

from PDB import ProteinStructure as PS
import copy
import numpy as np
import sys
import os


def refine_interactions(interaction_dict):
    '''
    :param interaction_dict: A dictionary with a identificator for both type of chains as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    :return: A dictionary with the same format of interaction_dict but with non-redundant interactions
    '''
    refined_interacciones = dict()
    for pair_type in interaction_dict:
        refined_interacciones.setdefault(pair_type, dict())
        for interaccion in interaction_dict[pair_type]:
            valid = True
            for otra_interaccion in refined_interacciones[pair_type]:
                intersect = set(interaccion).intersection(set(otra_interaccion))
                if len(intersect) >= max([len(interaccion), len(otra_interaccion)]) * 0.3:
                    if len(interaccion) >= len(otra_interaccion):
                        refined_interacciones.setdefault(pair_type, dict())[interaccion] = interaction_dict[pair_type][
                            interaccion]
                        del refined_interacciones[pair_type][otra_interaccion]
                        break
                    else:
                        valid = False
                        continue
            if interaccion not in refined_interacciones[pair_type] and valid:
                refined_interacciones.setdefault(pair_type, dict())[interaccion] = interaction_dict[pair_type][
                    interaccion]
    return refined_interacciones


def deconstruct_macrocomplex_by_interactions(pdb_file, output_folder='./', translation=np.array([0, 0, 0]),
                                             rotation=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])):
    '''
    :param pdb_file(str): A path to a pdb file
    :param output_folder (str): A path to a directory where the interactions will be stored
    :param translation (numpy array 3x1): a np array defining the translation of the new pdb's
    :param rotate (numpy array 3x3): a np array defining the rotation of the new pdb's
    '''
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    pdb = PS('original', pdb_file)
    pdb_list = list()
    todas_interacciones = dict()
    seqs = list()
    for chain in pdb:
        if chain.get_sequence_str() not in seqs:
            seqs.append(chain.get_sequence_str())
        for other_chain in pdb:
            if other_chain.get_sequence_str() not in seqs:
                seqs.append(other_chain.get_sequence_str())
            interacting_res_tuple = None
            if chain is not other_chain:
                int = chain.interacting_residues(other_chain, dist=3)
                if int:
                    int.extend(other_chain.interacting_residues(chain, dist=3))
                    interacting_res_tuple = tuple(int)
                if interacting_res_tuple:
                    todas_interacciones.setdefault(
                        ''.join(sorted(chain.get_sequence_str() + other_chain.get_sequence_str())), dict())[
                        interacting_res_tuple] = chain.get_id() + other_chain.get_id()
    interacciones_unicas = refine_interactions(todas_interacciones)
    pdb_list = list()
    for pair_type in interacciones_unicas:
        for interaccion in interacciones_unicas[pair_type]:
            pdb.transform(tran=translation, rot=rotation)
            pdb.save_to_file('%s/%s.pdb' % (output_folder, interacciones_unicas[pair_type][interaccion]),
                             chain_name=interacciones_unicas[pair_type][interaccion])
            pdb_list.append(PS(pdb_file='%s/%s.pdb' % (output_folder, interacciones_unicas[pair_type][interaccion]), id=interacciones_unicas[pair_type][interaccion]))
    return pdb_list
'''
Rotation matrix for a angle = A is the following:
        1       0       0
Rx(A) = 0       cosA    -sinA
        0       sinA    cosA
        cosA    0       sinA
Ry(A) = 0       1       0
        -sinA   0       cosA
        cosA    -sinA   0
Rz(A) = sinA    cosA    0
        0       0       1
                                                        '''

if __name__ == '__main__':

    # deconstruct_macrocomplex_by_interactions(sys.argv[1], sys.argv[2])
    deconstruct_macrocomplex_by_interactions(sys.argv[1], sys.argv[2])
