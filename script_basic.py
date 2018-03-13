#!/usr/bin/python3

from Bio.PDB import *
from sys import argv
import re

def superimpose(fix_chain,mov_chain):
	"""
	Superimpose peptide chain mov_chain on fix_chain with Biopython superimposer.
	"""
	sup = Superimposer()

	#//Obtenir una llista d'atoms-objects per a cada cadena
	fix_atoms = list(fix_chain.get_atoms())
	mov_atoms = list(mov_chain.get_atoms())
	#//Fer superimposicio
	sup.set_atoms(fix_atoms	,mov_atoms)

	#//Matriu amb Rotacio/translacio/rmsd a aplicar a la estructura
	print(sup.rotran)
	print(sup.rms)

	#//Aplicar rotacio/translacio als atoms mov. Per aplicar la mateixa rotacio a altres cadenes es podria fer servir la matriu anterior
	sup.apply(mov_atoms)

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

# Create a structure object. All structure objects are in a dictonary 
# We also obtain a set with the different chains (subunits of the macromolecula) 
# Moreover, we create a dictionary that contains a list with the id subunits interactiong by filename (key)

#Create a structure object for every file and store them as values in a dictionary with name of file as key
allpdb = { filename : parser.get_structure(filename, filename) for filename in argv if argv.index(filename) != 0}

#Store all chain id in a set. Create a dictionary where every filename has a list of chain-objects as values
subunits = set()
chainsbyfile = {}
for pdb in allpdb:
	chainsbyfile[pdb] = list()
	for subunit in allpdb[pdb].get_chains():
		chainsbyfile[pdb].append(subunit)
		subunits.add(subunit.get_id())

# Create a dictionary with all the subunits as keys and, as values, a list with the subunits which interacts with
###########NOTA D'OPIMITZACIO: Aixo es podria acoblar al loop anterior######################### 
interactions = {}
intlist = list()

for element in subunits: 
	for filename in allpdb:
		if element in chainsbyfile[filename]:
			interact = chainsbyfile[filename] 
			index = interact.index(element)
			if index == 0:
				intlist.append(interact[1])
				
			else:
				intlist.append(interact[0])
	interactions[element] = intlist
	intlist = list()

#
superimpose(chainsbyfile["1fn3_AB.pdb"][0], chainsbyfile['1fn3_AD.pdb'][0])