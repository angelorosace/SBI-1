#!/usr/bin/python3

from Bio.PDB import *
from sys import argv


p = PDBParser(PERMISSIVE = 1)

allpdb = {}
chains = {}
for filename in argv:
	if argv.index(filename) == 0:
		continue
	else:
		allpdb[filename] = p.get_structure(filename, filename)
		



for element in allpdb:
	for a in element:
		print(element["A"])
	#print(list(element.get_chains()))
