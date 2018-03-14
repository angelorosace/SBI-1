#!/usr/bin/python3

from Bio.PDB import *
from sys import argv
import re
import os

def superimpose(fix_chain,mov_chain,apply_chain):
	"""
	Superimpose peptide chain mov_chain on fix_chain with Biopython superimposer, and apply movement to chain .
	"""
	sup = Superimposer()

	#//Obtenir una llista d'atoms-objects per a cada cadena
	fix_atoms = list(fix_chain.get_atoms())
	mov_atoms = list(mov_chain.get_atoms())
	#//Fer superimposicio
	sup.set_atoms(fix_atoms	,mov_atoms)

	#//Matriu amb Rotacio/translacio/rmsd a aplicar a la estructura
	#print(sup.rotran)
	#print(sup.rms)

	sup.apply(apply_chain)
	return apply_chain

io = PDBIO()

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

##########################################
##Store Relations between files and chains
##########################################

"""
Dictionary and list index:
	1. allpdb:
		keys: input pdb filenames
		values: Structure pdb object from file

	2. Subunits: all chain names introduced
	3. Chainsbyfile: 
		keys: input pdb filenames
		values: list of chain-objects in pdb file
	4. Interactions: 
		key: chain name
		value: list of chain names wich interacts with

"""
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
###########NOTA A DAVID: Hi ha d'haver una forma millor de fer aixo######################### 
interactions = {}
intset = set()

for element in subunits: 
	for filename in allpdb:
		chainames = [ a.get_id() for a in chainsbyfile[filename] ]
		if element in chainames:
			index = chainames.index(element)
			if index == 0:
				intset.add(chainames[1])
				
			else:
				intset.add(chainames[0])
	interactions[element] = intset
	intset = set()

print(allpdb)
print(subunits)
print(chainsbyfile)
print(interactions)

"""
//NOTA A SANDVITX: He canviat els pdbs que utiltzem com a probes pels que fan servir Oriol i Alvaro. Perque? perque no hi havia caigut en que, malgrat estiguin en fitxers separats, les nostres
cadenes ja estan positionades on els toca, d'ons venen d'un fitxer comu
Ara surten un munt de warnings, pero no t'alarmis, funciona igualment. Es perque els pdbs que utilitzen Alvaro i Oriol els falta una columna 
"""

#//Sobre com fer el pdb output. No he trobat manera humana de fer-ho d'una forma mes neta. Si en trobes una millor, fes'la, pero no t'hi matis buscant. Ja he perdut jo prou temps com perque t'hi matis tu ara
io.set_structure(chainsbyfile["AB.pdb"][0])
io.save("temp1.pdb")
io.set_structure(chainsbyfile["AB.pdb"][1])
io.save("temp2.pdb")
moved_chain = superimpose(chainsbyfile["AB.pdb"][0], chainsbyfile['AD.pdb'][0], chainsbyfile['AD.pdb'][1])
io.set_structure(moved_chain)
io.save("temp3.pdb")
os.system("cat temp*.pdb > krosis.pdb")
