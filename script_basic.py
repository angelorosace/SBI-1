#!/usr/bin/python3

from Bio.PDB import *
from sys import argv


p = PDBParser(PERMISSIVE = 1)

# Create pdb objects. Save them in a dictionary. 
allpdb = {}
chains = set()

for filename in argv:
	if argv.index(filename) == 0:
		continue
	else:
		allpdb[filename] = p.get_structure(filename, filename)
		for chain in allpdb[filename].get_chains():
			chains.add(chain)
		#print(set((allpdb[filename].get_chains())).difference(chains))

		
		#chains.update(set((allpdb[filename].get_chains())).difference(chains))
		#print(chains)
		#chains.update(list(allpdb[filename].get_chains()))

# Obtain a dictionary with a list list of interactions (keys are chains)
#interactions = {}
#chains = []
#	for pdb in allpdb:
#		for chain in list(element.get_chains()):
#			interacions[chains] = append.chains 



#for element in allpdb:
#	for a in element:
#		print(element["A"])
	#print(list(element.get_chains()))

print(chains)
