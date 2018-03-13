#!/usr/bin/python3

from Bio.PDB import *
from sys import argv
import re

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

# Create a structure object. All structure objects are in a dictonary 
# We also obtain a set with the different chains (subunits of the macromolecula) 
# Moreover, we create a dictionary that contains a list with the id subunits interactiong by filename (key)
allpdb = {}
subunits = set()
chainsbyfile = {}
chainslist = list()

for filename in argv:
	if argv.index(filename) == 0:
		continue
	else:
		allpdb[filename] = parser.get_structure(filename, filename)
		subuid_last = str()
		for subunit in allpdb[filename].get_chains():
			suid = subunit.get_id()
			chainslist.append(suid)
			subunits.add(suid)
		chainsbyfile[filename] = chainslist
		chainslist = list()

			
# Create a dictionary with all the subunits as keys and a list with the subunits with which interacts
interactions = {}
intlist = list()

for element in subunits: 
	for filename in argv:
		if argv.index(filename) == 0:
			continue
		else:
			if element in chainsbyfile[filename]:
				interact = chainsbyfile[filename] 
				index = interact.index(element)
				if index == 0:
					intlist.append(interact[1])
					
				else:
					intlist.append(interact[0])
	interactions[element] = intlist
	intlist = list()



		


		
# Obtain a dictionary with a list of interactions (keys are subunits -chains in the case of simple complex-)
#interactions = {}
#chains = []
#for element in subunits:
#	for pdb in allpdb:
#		for chain in :
#			interacions[chains] = append.chains 


#for element in allpdb:
#	for a in element:
#		print(element["A"])
	#print(list(element.get_chains()))

#my_string =  "<Chain id=C>"
#p = re.compile("<Chain id=(.*)>")
#m = p.match(my_string)
#hey = m.groups()
#print(hey[0])

#print(chains)
#chains = set()
#for structure in allpdb:
#	#print(structure)
#	for model in allpdb[structure]:
#		#print(model)
#		for chain in model:
#			chains.add(chain.get_id())
			#print(dir(chain))
			#print(help(id))
			#m = p.match(chain)
			#chainid = m.groups()
			#chains.add(chainid[0])


#interactions[subuid] = intlist.append(subuid_last)
#			interactions[subuid_last] = intlist.append(subuid)
#			subuid_last = subuid
#		subuid_last = str()