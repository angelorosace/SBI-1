#!/usr/bin/python3

from Bio.PDB import *
from shutil import rmtree
from sys import exit
import os
import copy
import argparse

io = PDBIO()

#Counter of chains
global counter
counter = 1

def set_temp_name(number):
	"""
	Creates a name for a temporary file
	"""
	return str(options.path) + options.outfile + "_temp/temp" + str(number) + ".pdb"

def end_matchprot():
	"""
	Ends matchprot, but before it merges all temporary files into the final model and remove the temporary files
	"""
	#Merge temporary files into definitive result
	os.system(str("cat " + set_temp_name("*") + " > " + str(options.path) + str(options.outfile) + ".pdb"))

	#Delete temporary files (if not option temp is selected)
	if not options.temp:
		rmtree(str(options.path) + options.outfile + "_temp") 
	exit()


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
	#if clashcounter > 0:
	#	stderr.write("%i per cent atoms in chain %s are clashing\n" %(clashpercent, chain1.get_id()))

	#Return percentage of clashing atoms
	return clashpercent

def superimpose(fix_chain,mov_chain,apply_chain):
	"""
	Superimpose peptide chain mov_chain on fix_chain with Biopython superimposer, and return superimposition matrices
	"""
	sup = Superimposer()

	#Obtain atom-objects list by for fix and mov chain
	fix_atoms = list(fix_chain.get_atoms())
	mov_atoms = list(mov_chain.get_atoms())

	#Superimpose mov over fix
	sup.set_atoms(fix_atoms	,mov_atoms)

	#Do changes in applycain
	sup.apply(apply_chain)

	#Return superimposition matrices
	return sup

def recursive(chain_ob, camefrom_pdb):
	"""
	Main part of the script. 
	Checks all avalible interactions for chain_ob sequence. For each interaction, selects a chaintomove (the one that has the samesequence as chain_ob in the interaction)
	and superimposes it to chain_ob. The displacement matrices of this superimpositions are applied not to chaintomove, but to applychain (a deep copy ofthe chain that is
	interacting with chain to move). Applychain is displaced to interact with chain_ob.
	
	If applychain doesn't clash with any of the already displaced chain, it's saved on a temporary pdb file. Also, the recursive subrutine is re-executed with applychain
	as new chain_ob 
	
	If it clashes, is skipped and go to the next interacting chain of chain_ob
	"""

	#Extract id of chain object
	chain_ob_id = chain_ob.get_id()

	#Extract chain_ob ordinal (real sequence identifier, for )
	chain_ob_ordinal = idchainordinal[chain_ob_id]

	#Iterate over all the synonimes avalible for that chain
	for synonim in synonim_chains[chain_ob_ordinal]:

		#Iterate throught the pdbs where this sequene synonim is present
		for pdbinteracting in synonim_chains[chain_ob_ordinal][synonim]:

			#Iterate over chains in this pdbinteracting
			for chainteracting in chainsbyfile[pdbinteracting]:

				#Extract interacting chain id
				chainteracting_id = chainteracting.get_id()
				
				#Skip if interacting chain is same as synonim of chain_ob in the pdb
				if chainteracting_id == synonim:
					continue

				#Skip if interacting chain has already been placed
				if camefrom_pdb == pdbinteracting and chain_ob_id == synonim:
					continue

				#clash switch
				clashed_chain = False

				#Extract the pdb where chain_ob and chainteracting are both toghether
				pdbname = pdbinteracting
				pdb_ob = allpdb[pdbname]
				#chaintomove: chain of the same type than chain_ob in the pdb-file
				chaintomove = pdb_ob[0][synonim]
				#aplichain: chain that will be moved to interact with chain_ob
				applychain = copy.deepcopy(pdb_ob[0][chainteracting_id])

				#Verbose prints
				if options.verbose:
					print("Assembling chain %s from pdb %s to chain %s from assembled model." % (chainteracting_id, pdbname, chain_ob_id))

				#Superimpose
				supmatrix = superimpose(fix_chain = chain_ob, mov_chain = chaintomove, apply_chain = applychain)

				#Extract chain_applied ordinal (real sequence identifier, for )
				chain_applied_ordinal = idchainordinal[chainteracting_id]

				#Check if appliedchain collides with any of the already moved chains
				for chainmoved in coordsmoved[chain_applied_ordinal]:
					clashperc = clashtest(applychain, chainmoved)
					if clashperc > 80:
						clashed_chain = True
						#verbose print
						if options.verbose:
							print("Interacting chain ", chainteracting_id, " Skipped. Space already filled ")
						break

				#Check if applied chain has clashed with any of the already-moved chains. If so, reverse chain to original position and skip
				if clashed_chain == True:
					continue

				#Add one to chain limit if the ordinal of interacting chain corresponds to the limitant ordinal. Skip chain if limit is reached
				if (options.limitant_chain != "False") and (chain_applied_ordinal in limitant_ordinal.keys()):
					limitant_ordinal[chain_applied_ordinal] += 1
					if limitant_ordinal[chain_applied_ordinal] >= options.max_chains:
						continue

				#Add one to the superimpositions counter
				global counter
				counter += 1

				#Save appliedchain (the moved one) in a new pdb
				tempname =  set_temp_name(counter)
				savepdb(applychain, tempname)

				#Add applied atom coordenate list to coresponding list of already-moved chains
				coordsmoved[chain_applied_ordinal].append(get_chain_coords(applychain))

				#Verbose prints
				if options.verbose:
					print("Chain %s succesfully assembled to the model. Stored on temp %i" % (chainteracting_id,counter))

				#Check if limit of subunits has been reached, and end program if so (only when limitant chain option is not activated)
				if ( options.max_chains != -1) and (counter >= options.max_chains) and options.limitant_chain == "False":
					if options.verbose:
						print("maximum chains limit reached")
					end_matchprot()

				#Repeat process for appliedchain
				recursive(chain_ob = applychain, camefrom_pdb = pdbinteracting)

#########
##Options
#########

parser = argparse.ArgumentParser(description = "Matchprot reconstructs protein complexes from its individual protein interactions ")

parser.add_argument('-i', '--input',
	dest = "infiles",
	action = "store",
	nargs = '+',
	default = None, 
	help = "<Mandatory> Input PDB interaction files. Every file has to contain a unique, one-to-one, protein interaction")

parser.add_argument('-o', '--output',
	dest = "outfile",
	action = "store",    
	default = "protein_complex", 
	help = "(string) Output file name")

parser.add_argument('-d', '--path',
	dest = "path",
	action = "store",    
	default = "./", 
	help = "(string) Output directory")

parser.add_argument('-v', '--verbose',
	dest = "verbose",
	action = "store_true",    
	default = False, 
	help = "Print log comments in screen")

parser.add_argument('-t', '--temp',
	dest = "temp",
	action = "store_true",    
	default = False, 
	help = "Don't delete temporary files. Each one of them contain a subunit of the complex")

parser.add_argument('-m', '--max_chains',
	dest = "max_chains",
	action = "store",
	type = int,
	default = -1, 
	help = "(integer) Maximum number of subunits for the final model")

parser.add_argument('-l', '--limitant_chain',
	dest = "limitant_chain",
	action = "store",    
	default = "False", 
	help = "(string) If this option is activated, the --max_chains limit will only be applied to the subunits with the sequence of the selected identifier")


options = parser.parse_args()

#Print help and stop if there's no input 
if not options.infiles:
	print()
	parser.print_help()
	exit("Error: -i option is mandatory")

#Add "/" at the end of the path if user hasn't put it
if not options.path[-1] == "/":
	options.path = options.path + "/"

##########################################
##Store Relations between files and chains
##########################################

#//No es que no sapiga comptar, es que estic afegint i treient variables continuament i no valia la pena reindexar-ho tot
"""
Dictionary, list and other variables index:

	1. counter: counts the number of chains moved for bulding the present model

	1. allpdb:
		keys: input pdb filenames
		values: Structure pdb object from file

	2. chains_ids: all chain names introduced

	3. modelchains: list. Contains lists of atom coordenates for every already saved chain

	4. chainsbyfile: 
		keys: input pdb filenames
		values: list of chain-objects in pdb file

	5. originalseqs: dictionary
		keys: all different-sequence chain objects from the model
		values: "ordinal" identifier. Every unique sequence in the model has a unique ordinal identifier

	5. synonim_chains: 
		keys: all ordinals
			keys: chain identifiers with same ordinal (same sequence)
				values: pdb files names where this second chain is found

	7. idchainordinal: similar to original seqs. Dictionary: 
		keys: all chain identifiers in the model
		values: the corresponding ordinal sequence identifier

	7. coordsmoved: dictionary with
		keys: ordinals (unique sequence identifiers)
		values: chain_objects moved of the corresponding sequence type
"""

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

#Create allpdb dictionary
allpdb = { filename : parser.get_structure(filename, filename) for filename in options.infiles }

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

######################################################################
##Create synonim_chains, originalseqs and idchainordinal dictionaries
######################################################################

originalchains = dict()
synonim_chains = dict()
idchainordinal = dict()
countordinals = 1

#For pdb file
for pdbfile in chainsbyfile:

	#For testingchain by pdbfile, extract its seq and id
	for testingchain in chainsbyfile[pdbfile]:
		testingseq = ''.join(get_seq(testingchain))
		testingchainid = testingchain.get_id()

		#Compare with the original-chains dictionary, and classify the testingchainid with the corresponding ordinal
		for chain in originalchains:
			if comparechains(testingchain,chain):
				ordinal = originalchains[chain]

				if testingchainid in synonim_chains[ordinal].keys():
					synonim_chains[ordinal][testingchainid].append(pdbfile)
				else: 
					synonim_chains[ordinal][testingchainid] = [ pdbfile ]

				break

		#If this sequence hasn't a ordinal yet, create it
		else:
			ordinal = str(countordinals) + "th"
			originalchains[testingchain] = ordinal
			synonim_chains[ordinal] = dict()

			if testingchainid in synonim_chains[ordinal].keys():
				synonim_chains[ordinal][testingchainid].append(pdbfile)
			else: 
				synonim_chains[ordinal][testingchainid] = [ pdbfile ]

			countordinals += 1

		idchainordinal[testingchain.get_id()] = ordinal

#######################################################
##Setting coordsmoved and limitant ordinal dictionaries
#######################################################

#Create coordsmoved dictionary
coordsmoved = { chain: [] for chain in synonim_chains }

#Set limit chain dictionary control if option is activated
limitant_ordinal = dict()
if options.limitant_chain != "False":
	limitant_ordinal[idchainordinal[options.limitant_chain]] = 0 

#Debug prints
print("allpdb: ",allpdb)
print("chains_ids: ",chains_ids)
print("chainsbyfile: ",chainsbyfile)
print("coordsmoved: ",coordsmoved)
print("synonim_chains: ",synonim_chains)
print("original chains: ",originalchains)
print("idchainordinal :",idchainordinal)
print("limitant ordinal :",limitant_ordinal)

#####################################
##Initialize algorithm of matchmaking
#####################################

#Delete temporary files directory with same name if there was any
if os.path.isdir(str(options.path) + options.outfile + "_temp"):
	rmtree(str(options.path) + options.outfile + "_temp") 

#Make temporary files directory
os.mkdir(str(options.path) + options.outfile + "_temp")

#Obtain a random pdb interacting file
seed = list(allpdb.keys())[0]

#set starting variables
areclashes = False

#Iterate over pdb seed file chains
for chain in chainsbyfile[seed]:

	#Extract chain id
	chainid = chain.get_id()

	#Extract chain_ob ordinal (real sequence identifier, for )
	chain_ordinal = idchainordinal[chainid]

	for synonim in synonim_chains[chain_ordinal]:
		chainame = synonim
		#Check collisions for seed chains
		for chainmoved in coordsmoved[chain_ordinal]:
			clashperc = clashtest(chain, chainmoved)
			if clashperc > 80:
				areclashes = True 

		#save seed chains 
		if not areclashes:
			tempname =  set_temp_name(counter)
			savepdb(chain, tempname)
			coordsmoved[chain_ordinal].append(get_chain_coords(chain))

			#Start recursive
			recursive(chain_ob = chain, camefrom_pdb =  seed)
	prevchainid = chainid

#End program
end_matchprot()