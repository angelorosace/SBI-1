#!/usr/bin/python3

from Bio.PDB import *
from sys import argv
import os

def superimpose(fix_chain,mov_chain,apply_chain):
	"""
	Superimpose peptide chain mov_chain on fix_chain with Biopython superimposer, and apply movement to chain.
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
	1. Allpdb:
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
# The list is a list of dictionaries, becaise every chain id has a reference to the pdb were was found
###########NOTA A DAVID: Hi ha d'haver una forma millor de fer aixo######################### 
######## MARTA: He afegit d'on venen les chains (és a dir de quin pdb surten). Sinó era info que no et portava enlloc
interactions = {}
intset = list()
where = dict()

for element in subunits: 
	for filename in allpdb:
		chainames = [ a.get_id() for a in chainsbyfile[filename] ]
		if element in chainames:
			index = chainames.index(element)
			if index == 0:
				where[chainames[1]] = filename
				intset.append(where)
				
			else:
				where[chainames[0]] = filename
				intset.append(where)
				#intset.add(chainames[0])
			where = dict()
	interactions[element] = intset
	intset = list()


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
#//MARTA al habla: Ho deixo de mostra per no carregar-m'ho al intentar fer les iteracions per arribar a la macromolècula. 
#io.set_structure(chainsbyfile["AB.pdb"][0])
#io.save("temp1.pdb")
#io.set_structure(chainsbyfile["AB.pdb"][1])
#io.save("temp2.pdb")
#moved_chain = superimpose(chainsbyfile["AB.pdb"][0], chainsbyfile['AD.pdb'][0], chainsbyfile['AD.pdb'][1])
#io.set_structure(moved_chain)
#io.save("temp3.pdb")
#os.system("cat temp*.pdb > krosis.pdb")

# For each pdb we use the first chain as the fixed chain (reference) and save the second one
# We look for the other chains interacting with the fixed one and we move them
# Generate a pdb with the join interactions
# In each iteration we check if the superimposition has aleready been done, in that case we skip the file
predone = set()
done = set()
count = 2

for option in [[0,1],[1,0]]:
	for pdbname in chainsbyfile: 
		io.set_structure(chainsbyfile[pdbname][option[0]])
		io.save("temp1.pdb")
		io.set_structure(chainsbyfile[pdbname][option[1]])
		io.save("temp2.pdb")

		predone.add(pdbname)
		fixchain = chainsbyfile[pdbname][option[0]].get_id()
		interactwith = interactions[fixchain]
		for pdb2 in interactwith:
			print(pdb2)
			pdbinteract = list(pdb2.values())[0]
			print(pdbinteract)
			done = predone.copy()
			done.add(pdbinteract)
			if len(predone) == len(done): 
				continue
			else:
				count += 1
				testwhere = chainsbyfile[pdbinteract][option[0]].get_id()
				if testwhere == fixchain:
					posfix = option[0]
					posmov = option[1]
				else:
					posfix = option[1]
					posmov = option[0]
				print(pdbname, option[0], pdbinteract, posfix, pdbinteract, posmov)
				moved_chain = superimpose(chainsbyfile[pdbname][0], chainsbyfile[pdbinteract][posfix], chainsbyfile[pdbinteract][posmov])
				io.set_structure(moved_chain)
				tempname = "temp" + str(count) + ".pdb"
				io.save(tempname)
		predone = done
os.system("cat temp*.pdb > krosis.pdb")