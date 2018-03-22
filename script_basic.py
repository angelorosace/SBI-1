#!/usr/bin/python3

from Bio.PDB import *
from sys import *
import os
import copy

io = PDBIO()

#Counter of chains
global counter
counter = 0

def get_chain_coords(chain):
	"""
	Obtains atom coordenates list from a chain object
	"""
	return [ atom.get_coord() for atom in chain.get_atoms() ]

def get_seq(chain):
	"""
	Obtains sequence (list of AA identifiers) from a chain object
	"""
	return [ res.get_resname() for res in chain.get_residues() ]


def comparechains(chain1,chain2):
	"""
	Compares 2 sequences from 2 chain-objects and return true if they are equal (>95% matches)
	"""

	#Store AA sequence into lists
	listres1 = get_seq(chain1)
	listres2 = get_seq(chain2)

	#Compare if both lists are equal and store the number of matches and mismatches
	comparlist = [ res1 == res2 for res1,res2 in zip(listres1,listres2)]
	match_vs_mismatch = ((comparlist.count(True),comparlist.count(False)))

	#If there are more than 95% of matches, return true
	return (match_vs_mismatch[0]/len(comparlist) > 0.95)

def savepdb (pdb_ob,filename):
	"""
	Save an object in pdb format
	"""
	io.set_structure(pdb_ob)
	io.save(filename)

def clashtest(chain1, coordsatom2):
	"""
	Returns True if any of the atoms on chain1 clash with any of the atom coordenates in chain2
	"""
	clashcounter = 0

	#Obtain a list of atom objects for chain1
	atom1 = [ atom for atom in chain1.get_atoms() ]

	#Obtain number of atoms in both chains
	numatoms1 = len(atom1)
	numatoms2 = len(coordsatom2)

	#Prepare list of atoms for neighbour search
	ns = NeighborSearch(atom1)

	#Check collisions of every atom coordenates in coordsatom2 to already-prepared-list of atoms1 (collision distance: 4 amstrong)
	for clashes in map(lambda x: ns.search(x,4.0), coordsatom2):
		#Count number of atoms clashing
		if len(clashes)> 0:
			clashcounter += 1

	#If any atom clashes, calculate the overall percentage of clashing atoms in chain
	#//El missatge d'error el deixo de moment com a demostracio. A la versio definitiva no hi hauria de ser
	clashpercent  = clashcounter/numatoms1 * 100
	if clashcounter > 0:
		stderr.write("%i per cent atoms in chain %s are clashing\n" %(clashpercent, chain1.get_id()))

	#Return percentage of clashing atoms
	return clashpercent

def superimpose(fix_chain,mov_chain,apply_chain):
	"""
	Superimpose peptide chain mov_chain on fix_chain with Biopython superimposer, and apply movement to chain.
	"""
	sup = Superimposer()

	#Obtain atom-objects list by for fix and mov chain
	fix_atoms = list(fix_chain.get_atoms())
	mov_atoms = list(mov_chain.get_atoms())

	#Superimpose mov over fix
	sup.set_atoms(fix_atoms	,mov_atoms)

	#Apply superimposition movement to apply chain, and return it
	sup.apply(apply_chain)
	return apply_chain

def recursive(chain_ob, camefrom, interactions, allpdb):
	"""
	NOTA A MARTA: Aquesta subrutina, benvolguda Marta, potser ens solucioni gran part del treball o potser ens peti l'ordinador. L'important es que sortirem de dubtes
	L'he comentat molt perque la puguis entendre facil i rapidament
	"""

	#//Tancament de seguretat
	#if counter > 3:
	#	exit()

	#Extract chain_ob numeral (real sequence identifier, for )
	for originalchain in originalchains:
		if comparechains(originalchain, chain_ob):
			chain_ob_numeral = originalchains[originalchain]

	#Iterate over all the synonimes avalible for that chain
	for synonim in synonim_chains[chain_ob_numeral]:
		chainame = synonim
		#Iterate throught the chains which interact with the present chain
		for chainteracting in interactions[chainame]:

			#for every pdb where there's an interaction between chainteracting and chain_ob
			for pdbinteracting in interactions[chainame][chainteracting]:

				#clash switch
				clashed_chain = False

				#If this loop corresponds to the already superimposed pdb for this chain where the present chain_ob comes from, skip it
				if pdbinteracting == camefrom:
					continue

				#Extract the pdb where chain_ob and chainteracting are both toghether
				pdbname = pdbinteracting
				pdb_ob = allpdb[pdbname]
				#chaintomove: chain of the same type than chain_ob in the pdb-file
				chaintomove = pdb_ob[0][chainame]
				#aplichain: chain that will be moved to interact with chain_ob
				applychain = copy.deepcopy(pdb_ob[0][chainteracting])

				#Superimpose
				appliedchain = superimpose(fix_chain = chain_ob, mov_chain = chaintomove, apply_chain = applychain)

				#Extract chain_applied numeral (real sequence identifier, for )
				for originalchain in originalchains:
					if comparechains(originalchain, applychain):
						chain_applied_numeral = originalchains[originalchain]

				#Check if appliedchain collides with any of the already moved chains
				for chainmoved in coordsmoved[chain_applied_numeral]:
					clashperc = clashtest(appliedchain, chainmoved)
					if clashperc > 80:
						clashed_chain = True
						break

				#Check if applied chain has clashed with any of the already-moved chains, and skip chain  and pdb if so
				if clashed_chain == True:
					continue

				#Add one to the superimpositions counter
				global counter
				counter += 1

				#Save appliedchain (the moved one) in a new pdb
				tempname = "temp" + str(counter) + ".pdb"
				savepdb(appliedchain, tempname)

				#Add applied chain to coresponding list of already-moved chains
				coordsmoved[chain_applied_numeral].append(get_chain_coords(appliedchain))

				#Repeat process for appliedchain
				recursive(appliedchain, pdbinteracting, interactions, allpdb)


##########################################
##Store Relations between files and chains
##########################################

"""
Dictionary, list and other variables index:
	6. counter: counts the number of chains moved for bulding the present model
	1. Allpdb:
		keys: input pdb filenames
		values: Structure pdb object from file

	2. chains_ids: all chain names introduced
	3. modelchains: list. Contains lists of atom coordenates for every already saved chain
	4. Chainsbyfile: 
		keys: input pdb filenames
		values: list of chain-objects in pdb file
	5. Original chains: dictionary
		keys: all different chain objects from the model
		values: "numeral" identifier. Every unique sequence in the model has a unique numeral identifier
	5. synonim_chains: 
		keys: all numerals
			keys: chain identifiers with same numeral (same sequence)
				values: pdb files names where this second chain is found 
	5. Interactions: 
		key: chain name
		value: dictionary
			key: chains which key interact with
			value: pdbfile where this interaction is found
	7. coordsmoved: dictionary with
		keys: Numerals (unique sequence identifiers)
		values: chain_objects moved of the corresponding sequence type
"""

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

#Create allpdb dictionary
allpdb = { filename : parser.get_structure(filename, filename) for filename in argv if argv.index(filename) != 0}

#Create chains, chains_ids set and chainsbyfile dictionary
chains_ids = set()
modelchains = set()
chainsbyfile = {}
for pdb in allpdb:
	chainsbyfile[pdb] = set()
	for subunit in allpdb[pdb].get_chains():
		chains_ids.add(subunit.get_id())
		modelchains.add(subunit)
		chainsbyfile[pdb].add(subunit)

##################################################
##Create synonim_chains and originalseqs dictionary
##################################################

originalchains = dict()
synonim_chains = dict()
countnumerals = 1

#For pdb file
for pdbfile in chainsbyfile:

	#For testingchain by pdbfile, extract its seq and id
	for testingchain in chainsbyfile[pdbfile]:
		testingseq = ''.join(get_seq(testingchain))
		testingchainid = testingchain.get_id()

		#Compare with the original-chains dictionary, and classify the testingchainid with the corresponding numeral
		for chain in originalchains:
			if comparechains(testingchain,chain):
				numeral = originalchains[chain]
				synonim_chains[numeral][testingchainid] = pdbfile
				break

		#If this sequence hasn't a numeral yet, create it
		else:
			numeral = str(countnumerals)+"th"
			originalchains[testingchain] = numeral
			synonim_chains[numeral] = dict()
			synonim_chains[numeral][testingchainid] = pdbfile
			countnumerals += 1

################################
#Create intersections dictionary
################################

interactions = {}
intdict = dict()

#Iterate over every subunit in the model
for element in chains_ids:

	#Iterate over every file in input
	for filename in allpdb:
		chainames = [ a.get_id() for a in chainsbyfile[filename] ]

		#If this file contains the subunit of interest 
		if element in chainames:
			index = chainames.index(element)

			#Iterate over chains_ids in pdbfile
			for intchain in chainames:

				#Skip key chain
				if intchain == element:
					continue

				#Store interactions
				if intchain not in intdict.keys():
					intdict[intchain] = [ filename ]
				else:
					intdict[intchain].append(filename)
					
	interactions[element] = intdict
	intdict = dict()

#Create coordsmoved dictionary
coordsmoved = { chain: [] for chain in synonim_chains }

#print("allpdb: ",allpdb)
#print("chains_ids: ",chains_ids)
#print("chainsbyfile: ",chainsbyfile)
#print("interactions: ",interactions)
print("coordsmoved: ",coordsmoved)
print("synonim_chains: ",synonim_chains)
print("original chains: ",originalchains)

#####################################
##Initialize algorithm of matchmaking
#####################################

#Obtain a random pdb interacting file
seed = list(allpdb.keys())[0]
#Iterate over pdb seed file chains
areclashes = False
for chain in chainsbyfile[seed]:

	#Extract chain_ob numeral (real sequence identifier, for )
	for originalchain in originalchains:
		if comparechains(originalchain, chain):
			chain_numeral = originalchains[originalchain]

	#Check collisions for seed chains
	for chainmoved in coordsmoved[chain_numeral]:
		clashperc = clashtest(chain, chainmoved)
		if clashperc > 80:
			areclashes = True 

	#save seed chains 
	if not areclashes:
		savepdb(chain, "temp" + str(counter) + ".pdb")
		coordsmoved[chain_numeral].append(get_chain_coords(chain))

		recursive(chain, seed, interactions, allpdb)

os.system("cat temp*.pdb > krosis.pdb")

#Delete temporary files
os.system("rm temp*.pdb")
