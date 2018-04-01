#!/usr/bin/python3

from Bio.PDB import *
from shutil import rmtree
from sys import stdout
from os import system,mkdir
from os.path import isdir
import argparse

io = PDBIO()

#Counter of chains
global counter
counter = 1

def set_temp_name(number):
	"""
	Creates a name for a temporary file based on counter
	"""
	return str(options.path) + options.outfile + "_temp/temp" + str(number) + ".pdb"

def end_script():
	"""
	Ends script, but before it merges all temporary files into the final model and remove the temporary files
	"""
	#Merge temporary files into definitive result
	system(str("cat " + set_temp_name("*") + " > " + str(options.path) + str(options.outfile) + ".pdb"))

	#Delete temporary files (if not option temp is selected)
	if not options.temp:
		rmtree(str(options.path) + options.outfile + "_temp")
	exit()

def get_chain_coords_CA_P(chain):
	"""
	Returns a list of coordenates from the main atoms (C alpha for protein and P for DNA/RNA) of a chain
	"""
	return [ atom.get_coord() for atom in chain.get_atoms() if atom.get_id() == "CA" or atom.get_id() == "P"]

def get_seq(chain):
	"""
	Obtains sequence (list of AA identifiers) from a chain object
	"""
	return [ res.get_resname() for res in chain.get_residues() ]


def comparechains(chain1,chain2):
	"""
	Compares 2 sequences from 2 chain-objects and return True if they are equal (>95% matches)
	"""

	#Store AA sequence into lists
	listres1 = get_seq(chain1)
	listres2 = get_seq(chain2)

	#If length is not the same, sequences are considered not equal
	if len(listres1) != len(listres2):
		return False

	#Compare if both lists are equal and store the number of matches and mismatches
	comparlist = [ res1 == res2 for res1,res2 in zip(listres1,listres2)]
	match_vs_mismatch = ((comparlist.count(True),comparlist.count(False)))

	#If there are more than 95% of matches, return true
	return (match_vs_mismatch[0]/len(comparlist) > 0.95)

def savepdb (pdb_ob,filename):
	"""
	Save an object (mostly a chain) in pdb format
	"""
	io.set_structure(pdb_ob)
	io.save(filename)

def clashtest(atoms1, coordsatom2):
	"""
	Returns True if more than 10% of atoms in atoms1 are clashing with any atom in the coordenates of coordsatom2
	"""
	clashcounter = 0

	#Obtain number of atoms in both chains
	numatoms1 = len(atoms1)

	#Prepare list of atoms for neighbour search
	ns = NeighborSearch(atoms1)

	#Check collisions of every atom coordenates in coordsatom2 to already-prepared-list of atoms1 (collision distance: 4 amstrong)
	for clashes in map(lambda x: ns.search(x,2.0), coordsatom2):
		#Count number of atoms clashing
		if len(clashes)> 0:
			clashcounter += 1

	#If any atom clashes, calculate the overall percentage of clashing atoms in chain
	clashpercent  = clashcounter/numatoms1 * 100
	
	#Return True or False depending if more than 5% of atoms are clashing
	return clashpercent > 5

def superimpose(fix_chain,mov_chain,apply_chain):
	"""
	Superimpose mov_chain on fix_chain with Biopython superimposer, and apply the rotation and movement matrices to apply_chain
	"""

	sup = Superimposer()

	#Obtain list of main atom-objects (Carbon alpha for protein and Phosphate for DNA/RNA) by for fix and mov chain
	fix_atoms = list([atom for atom in fix_chain.get_atoms() if atom.get_id() == "CA" or atom.get_id() == "P"])
	mov_atoms = list([atom for atom in mov_chain.get_atoms() if atom.get_id() == "CA" or atom.get_id() == "P"])

	#Superimpose mov over fix
	sup.set_atoms(fix_atoms	,mov_atoms)

	#Do changes in applycain
	sup.apply(apply_chain)

	#Return superimposition matrices
	return sup

def recursive(chain_ob, camefrom_pdb):
	"""
	Main part of the script. 
	Checks all avalible interactions for chain_ob sequence. For each interaction, selects a chaintomove 
	(the one that has the samesequence as chain_ob in the interaction) and superimposes it to chain_ob.
	The displacement matrices of this superimpositions are applied not to chaintomove, but to applychain
	(a deep copy of the chain that is interacting with chain to move). Applychain is displaced to interact with chain_ob.
	
	If applychain doesn't clash with any of the already displaced chains in coordsmoved, it's saved on a
	temporary pdb file. Also, the recursive subrutine is re-executed with applychain as new chain_ob 
	
	If it clashes, applychain is skipped without saving and go to the next interacting chain of chain_ob
	"""

	#Extract id of chain object
	chain_ob_id = chain_ob.get_id()

	#Extract chain_ob ordinal
	chain_ob_ordinal = idchainordinal[chain_ob_id]

	#Iterate over all the synonimes avalible for that chain (different input ids for this sequence)
	for synonim in synonim_chains[chain_ob_ordinal]:

		#if unique option is activated, only use pdbs with the original chain id of this chain
		if options.unique and (synonim != chain_ob_id):
			continue

		#Iterate throught the pdbs where this sequence synonim is present
		for pdbinteracting in synonim_chains[chain_ob_ordinal][synonim]:

			#Iterate throught the models  of the pdb where synonim is present
			for model_syn in synonim_chains[chain_ob_ordinal][synonim][pdbinteracting]:

				#Iterate over chains in this model
				for model_inter in chainsbyfilebymodel[pdbinteracting]:

					for chainteracting in chainsbyfilebymodel[pdbinteracting][model_inter]:

						#Extract interacting chain id
						chainteracting_id = chainteracting.get_id()
						
						#Skip if interacting chain and chain to move are actually the same
						if (synonim == chainteracting_id) and (model_syn == model_inter):
							continue

						#Skip if interacting chain has already been placed in previous interaction
						if (camefrom_pdb == pdbinteracting) and (chain_ob_id == chainteracting):
							continue

						#clash switch
						clashed_chain = False

						#Extract the pdb where chain_ob and chainteracting are both toghether
						pdb_ob = allpdb[pdbinteracting]
						#chaintomove: chain of the same type than chain_ob in the pdb-file-model 
						chaintomove = pdb_ob[model_syn][synonim]
						#aplichain: chain that will be moved to interact with chain_ob
						applychain = pdb_ob[model_inter][chainteracting_id].copy()

						#Verbose prints
						if options.verbose:
							stdout.write("Assembling chain %s from pdb %s to chain %s from assembled model.\n" % (chainteracting_id, pdbinteracting, chain_ob_id))

						#Superimpose
						supmatrix = superimpose(fix_chain = chain_ob, mov_chain = chaintomove, apply_chain = applychain)

						#Extract chain_applied ordinal (real sequence identifier, for )
						chain_applied_ordinal = idchainordinal[chainteracting_id]

						#Check if appliedchain collides with any of the already moved chains
						applyatoms = [ atom for atom in applychain.get_atoms() if atom.get_id() == "CA" or atom.get_id() == "P"]
						for counter in coordsmoved:
							if clashtest(applyatoms, coordsmoved[counter]):
								clashed_chain = True
								#verbose print
								if options.verbose:
									stdout.write("Interacting chain " + str(chainteracting_id) + " Skipped. Space already filled \n") 
								break

						#Check if applied chain has clashed with any of the already-moved chains. If so, reverse chain to original position and skip
						if clashed_chain == True:
							continue

						#If unique option is activated use chain id as keys, else use ordinals-sequence
						if options.unique:
							mykey = chainteracting_id
						else: 
							mykey = chain_applied_ordinal

						#If limitant chains is activated, add one to chain limit if the ordinal of interacting chain corresponds to the limitant ordinal. Skip chain if limit is reached
						if (options.limitant_chains != "False") and (mykey in limitant_ordinals.keys()):
							limitant_ordinals[mykey][0] += 1
							if limitant_ordinals[mykey][0] > limitant_ordinals[mykey][1]:
								continue

						#Add one to the superimpositions counter
						counter = counter + 1

						#Save appliedchain (the moved one) in a new pdb with a tempname
						tempname =  set_temp_name(counter)
						savepdb(applychain, tempname)

						#Add main-atom coordinates from applychain to coordsmoved 
						coordsmoved[counter] = get_chain_coords_CA_P(applychain)

						#Verbose prints
						if options.verbose:
							stdout.write("Chain %s succesfully assembled to the model. Stored on temp %i\n" % (chainteracting_id,counter))

						#Check if limit of subunits has been reached, and end program if so (only when limitant chain option is not activated)
						if ( options.max_chains != -1) and (counter >= options.max_chains[0]) and options.limitant_chains == "False":
							if options.verbose:
								stdout.write("maximum chains limit reached")
							end_script()

						#Repeat process for appliedchain
						recursive(chain_ob = applychain, camefrom_pdb = pdbinteracting)

#########
##Options
#########

parser = argparse.ArgumentParser(description = "script reconstructs protein complexes from its individual protein interactions ")

parser.add_argument('-i', '--input',
	dest = "infiles",
	action = "store",
	nargs = '+',
	default = None, 
	help = "<Mandatory> Input PDB interaction files. Do not use the same id for different sequences")

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
	help = "print log comments in stdout")

parser.add_argument('-t', '--temp',
	dest = "temp",
	action = "store_true",    
	default = False, 
	help = "Don't delete temporary files. Each one of them contain a subunit of the output complex")

parser.add_argument('-m', '--max_chains',
	dest = "max_chains",
	action = "store",
	type = int,
	default =  -1,
	nargs = '+', 
	help = "(integer) Maximum number of subunits for the final model. If limitant_chains is activated, absolute sti")

parser.add_argument('-l', '--limitant_chains',
	dest = "limitant_chains",
	action = "store",    
	default = "False",
	nargs = '+', 
	help = """(string) List of limitant chain identifiers for the model.
	Each sequence specified in this list will be in the model equal or less times than the corresponding number in --limitant_chains list multiplied by the --proportions_multiplier option.
	If a subunit has more than one name in your model, put just one of its identifiers.""")

parser.add_argument('-p', '--proportions_multiplier',
	dest = "proportions_multiplier",
	action = "store",
	type = int,
	default = 1,
	help = "(integer) If --limitant_chains is activated, number to multiply the proportions by. Default is 1")

parser.add_argument('-s', '--seed_pdb',
	dest = "seed_pdb",
	action = "store",
	default = "False",
	help = "(string) Name of the pdb input file from which the reconstruction will start. Default is random")

parser.add_argument('-u', '--unique_ids',
	dest = "unique",
	action = "store_true",
	default = False,
	help = """If this option is activated, all input chains with different chain ids will be treated as different subunits,
	 regardless of the sequence composition""")

options = parser.parse_args()

#Print help and stop if there's no input 
if not options.infiles:
	print()
	parser.print_help()
	exit("Error: -i option is mandatory")

#If limitant_chains option is activated and its length do not coincide with limitant_chains, exit and show error message
if (options.limitant_chains != "False") and (len(options.limitant_chains) != len(options.max_chains)):
	print()
	parser.print_help()
	exit("Error: max_chains list of integers and limitant_chainss list of chains are not the same length")

#Add "/" at the end of the path if user hasn't put it
if not options.path[-1] == "/":
	options.path = options.path + "/"

##########################################
##Store Relations between files and chains
##########################################

"""
Dictionary, list and other variables index:

	1. counter: counts the number of chains moved for bulding the present model

	2. allpdb:
		keys: input pdb filenames
		values: Structure pdb object from file

	3. chainsbyfilebymodel: 
		keys: input pdb filenames
			keys: model ids of each pdbfile
			values: list of chain-objects in model of pdb file

	4. originalseqs: dictionary
		keys: all different-sequence chain objects from the model
		values: "ordinal" identifier. Every unique sequence in the model has a unique ordinal identifier

	5. synonim_chains: 
		keys: all ordinals
			keys: chain identifiers with same ordinal (same sequence)
				keys: pdb files names where this chain is found
					values: list of models from the pdbfiles where chains are found

	6. idchainordinal: similar to original seqs. Dictionary: 
		keys: all chain identifiers in the model
		values: the corresponding ordinal sequence identifier

	7. coordsmoved: set with the Calpha or P coordinates of all chains already saved 

	8. limitant_ordinals: dictinoary. Empty if limitant chains option is not activated
		keys: all ordinals for chains stated in limitant_chains option
			value[0]: counter of how many sequences corresponding to this ordinal have been already saved
			value[1]: maximum number allowed of chains with this sequence in the final model
"""

#Verbose prints
if options.verbose:
	stdout.write("Processing input...\n")

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

#Create allpdb dictionary
allpdb = { filename : parser.get_structure(filename, filename) for filename in options.infiles }

######################################
#Create chainsbyfilebymodel dictionary
######################################

chainsbyfilebymodel = {}

#For pdb file
for pdb in allpdb:

	chainsbyfilebymodel[pdb] = dict()

	#For model in pdb file
	for model in allpdb[pdb].get_models():
		model_id = model.get_id()
		chainsbyfilebymodel[pdb][model_id] = set()

		#For chain in model of pdb file
		for subunit in allpdb[pdb][model_id].get_chains():
			chainsbyfilebymodel[pdb][model_id].add(subunit)

######################################################################
##Create synonim_chains, originalseqs and idchainordinal dictionaries
######################################################################

originalseqs = dict()
synonim_chains = dict()
idchainordinal = dict()
countordinals = 1

#For pdb file
for pdbfile in chainsbyfilebymodel:

	#For model in pdb file
	for model in chainsbyfilebymodel[pdbfile]:

		#For testingchain by model, extract its seq and id
		for testingchain in chainsbyfilebymodel[pdbfile][model]:
			testingseq = ''.join(get_seq(testingchain))
			testingchainid = testingchain.get_id()

			#Compare with the original-chains dictionary, and classify the testingchainid with the corresponding ordinal
			for chain in originalseqs:
				if comparechains(testingchain,chain):
					ordinal = originalseqs[chain]

					#If its a new synonim
					if testingchainid in synonim_chains[ordinal].keys():

						#if model for this synonim is already seen 
						if pdbfile in synonim_chains[ordinal][testingchainid].keys():
							synonim_chains[ordinal][testingchainid][pdbfile].append(model)

						else:
							synonim_chains[ordinal][testingchainid][pdbfile] = [ model ]

					#If we see this synonim for first time
					else: 
						synonim_chains[ordinal][testingchainid] = { pdbfile : [ model ] } 

					break

			#If this sequence hasn't a ordinal yet, create it and add present chain, pdb, model and so
			else:
				ordinal = str(countordinals) + "th"
				originalseqs[testingchain] = ordinal
				synonim_chains[ordinal] = dict()

				if testingchainid in synonim_chains[ordinal].keys():
					synonim_chains[ordinal][testingchainid][pdbfile].append(model)
				else: 
					synonim_chains[ordinal][testingchainid] = { pdbfile : [ model ] }

				countordinals += 1

			idchainordinal[testingchain.get_id()] = ordinal

#######################################################
##Setting coordsmoved and limitant ordinal dictionaries
#######################################################

#Create coordsmoved dictionary
coordsmoved = {}

#Set limitant ordinals double dictionary if limitant_chains option is set
limitant_ordinals = dict()
if options.limitant_chains != "False":

	#if unique option names is activated, use chane ids as keys instead of sequence-ordinals
	if options.unique:
		limitant_ordinals = { chainid : [ 0, chainlimit*options.proportions_multiplier ] for chainid,chainlimit in zip(options.limitant_chains,options.max_chains) }

	#Else, use ordinals as keys
	else: 
		limitant_ordinals = { idchainordinal[chainid] : [ 0, chainlimit*options.proportions_multiplier ] for chainid,chainlimit in zip(options.limitant_chains,options.max_chains) }

#verbose proportions_multiplier
if options.verbose:
	stdout.write("\nnumber of input pdbs: " + str(len(allpdb)) + "\n")
	stdout.write("number of input chains (with different id): " + str(len(idchainordinal)) + "\n")
	stdout.write("number of unique chain sequences: " + str(len(originalseqs)) + "\n\n")

#Debug prints
"""
print("allpdb: ",allpdb)
print("chainsbyfilebymodel: ",chainsbyfilebymodel)
print("coordsmoved: ",coordsmoved)
print("synonim_chains: ",synonim_chains)
print("original chains: ",originalseqs)
print("idchainordinal :",idchainordinal)
print("limitant ordinals :",limitant_ordinals)
"""

#####################################
##Initialize algorithm of matchmaking
#####################################

#Delete temporary files directory with same name if there was any
if isdir(str(options.path) + options.outfile + "_temp"):
	rmtree(str(options.path) + options.outfile + "_temp") 

#Make temporary files directory
mkdir(str(options.path) + options.outfile + "_temp")

#Obtain a random pdb interacting file
if options.seed_pdb != "False":
	seed = options.seed_pdb
else:
	seed = list(allpdb.keys())[0]

#set starting variables
areclashes = False

#Iterate over pdb seed file models
for model in chainsbyfilebymodel[seed]:

	#Iterate over chains of the model
	for chain in chainsbyfilebymodel[seed][model]:

		#Extract chain id
		chainid = chain.get_id()

		#Extract chain_ob ordinal (real sequence identifier, for )
		chain_ordinal = idchainordinal[chainid]

		for synonim in synonim_chains[chain_ordinal]:
			chainame = synonim

			#if unique option is activated, only use pdbs with the original chain id of this chain
			if options.unique and (synonim != chainid):
				continue

			#If unique option is activated use chain is as keys, else use ordinals-sequence
			if options.unique:
				mykey = chainid
			else: 
				mykey = chain_ordinal

			#Add one to chain limit if the ordinal of interacting chain corresponds to the limitant ordinal. Skip chain if limit is reached
			if (options.limitant_chains != "False") and (mykey in limitant_ordinals.keys()):
				limitant_ordinals[mykey][0] += 1
				if limitant_ordinals[mykey][0] > limitant_ordinals[mykey][1]:
					continue

			#Check if limit of subunits has been reached, and end program if so (only when limitant chain option is not activated)
			if ( options.max_chains != -1) and (counter >= options.max_chains[0]) and options.limitant_chains == "False":
				if options.verbose:
					stdout.write("maximum chains limit reached")
				end_script()

			#Check collisions for seed chains
			chainatoms = [ atom for atom in chain.get_atoms() if atom.get_id() == "CA" or atom.get_id() == "P"]
			for counter in coordsmoved:
				if clashtest(chainatoms, coordsmoved[counter]):
					areclashes = True 

			#save seed chains
			if not areclashes:
				tempname =  set_temp_name(counter)
				savepdb(chain, tempname)
				coordsmoved[counter] = (get_chain_coords_CA_P(chain))

				#Start recursive
				recursive(chain_ob = chain, camefrom_pdb =  seed)

#End program
end_script()
