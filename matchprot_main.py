#!/usr/bin/python3

from Bio.PDB import *
from shutil import rmtree
from sys import *
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
	Save an object in pdb format
	"""
	io.set_structure(pdb_ob)
	io.save(filename)

def clashtest(chain1, coordsatom2):
	"""
	Returns True if more than 10% of atoms in chain1 are clashing with any atom in the coordenates of coordsatom2
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
	return clashpercent > 10

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
				for counter in coordsmoved:
					if clashtest(applychain, coordsmoved[counter]):
						clashed_chain = True
						#verbose print
						if options.verbose:
							print("Interacting chain ", chainteracting_id, " Skipped. Space already filled ")
						break

				#Check if applied chain has clashed with any of the already-moved chains. If so, reverse chain to original position and skip
				if clashed_chain == True:
					continue

				#Add one to chain limit if the ordinal of interacting chain corresponds to the limitant ordinal. Skip chain if limit is reached
				if (options.limitant_chains != "False") and (chain_applied_ordinal in limitant_ordinals.keys()):
					limitant_ordinals[chain_applied_ordinal][0] += 1
					if limitant_ordinals[chain_applied_ordinal][0] > limitant_ordinals[chain_applied_ordinal][1]:
						continue

				#Add one to the superimpositions counter
				counter = counter + 1

				#Save appliedchain (the moved one) in a new pdb
				tempname =  set_temp_name(counter)
				savepdb(applychain, tempname)

				#Add applied atom coordenate list to coresponding list of already-moved chains
				coordsmoved[counter] = get_chain_coords(applychain)

				#Verbose prints
				if options.verbose:
					print("Chain %s succesfully assembled to the model. Stored on temp %i" % (chainteracting_id,counter))

				#Check if limit of subunits has been reached, and end program if so (only when limitant chain option is not activated)
				if ( options.max_chains != -1) and (counter >= options.max_chains[0]) and options.limitant_chains == "False":
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
	default =  -1,
	nargs = '+', 
	help = "(integer) Maximum number of subunits for the final model. If limitant_chains is activated, list of absolute proportions for subunit in complex")

parser.add_argument('-l', '--limitant_chains',
	dest = "limitant_chains",
	action = "store",    
	default = "False",
	nargs = '+', 
	help = """(string) List of limitant chain identifiers for the model.
	Each sequence specified in this list will be in the model equal or less times than the corresponding number in --limitant_chains list multiplied by the --proportions_multiplier option.
	If a chain has more than one name in your model, put just one of its identifiers.""")

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
	help = "(string) Name of the pdb file from which the reconstruction will start. Default is random")

options = parser.parse_args()

#Print help and stop if there's no input 
if not options.infiles:
	print()
	parser.print_help()
	exit("Error: -i option is mandatory")

#If limitant_chains option is activated and its length do not coincide with limitant_chains, die and show error message
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

	7. coordsmoved: set with the coordinates of all chains already saved 

	8. limitant_ordinals: dictinoary. Empty if limitant chains option is not activated
		keys: all ordinals for chains stated in limitant_chains option
			value[0]: counter of how many sequences corresponding to this ordinal have been already saved
			value[1]: maximum number allowed of chains with this sequence in the final model
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
coordsmoved = {}

#Set limitant ordinals double dictionary if limitant_chains option is set
limitant_ordinals = dict()
if options.limitant_chains != "False":
	limitant_ordinals = { idchainordinal[chainid] : [ 0, chainlimit*options.proportions_multiplier ] for chainid,chainlimit in zip(options.limitant_chains,options.max_chains) }

#Debug prints
print("allpdb: ",allpdb)
print("chains_ids: ",chains_ids)
print("chainsbyfile: ",chainsbyfile)
print("coordsmoved: ",coordsmoved)
print("synonim_chains: ",synonim_chains)
print("original chains: ",originalchains)
print("idchainordinal :",idchainordinal)
print("limitant ordinals :",limitant_ordinals)

#####################################
##Initialize algorithm of matchmaking
#####################################

#Delete temporary files directory with same name if there was any
if os.path.isdir(str(options.path) + options.outfile + "_temp"):
	rmtree(str(options.path) + options.outfile + "_temp") 

#Make temporary files directory
os.mkdir(str(options.path) + options.outfile + "_temp")

#Obtain a random pdb interacting file
if options.seed_pdb != "False":
	seed = options.seed_pdb
else:
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
		for counter in coordsmoved:
			if clashtest(chain, coordsmoved[counter]):
				areclashes = True 

		#save seed chains 
		if not areclashes:
			tempname =  set_temp_name(counter)
			savepdb(chain, tempname)
			coordsmoved[counter] = (get_chain_coords(chain))

			#Start recursive
			recursive(chain_ob = chain, camefrom_pdb =  seed)
	prevchainid = chainid

#End program
end_matchprot()

#########################################################
##########HOLA MARTAAAAAAAAAAAAAAAAAAAAAAAA: Log change
#######################################################

#//Per fer front al repte d'ensamblar una estructura asimètrica (com es l'ATP syntasa), he fet les seguents modificacions

#//1. Clashtest global: Ara el clash test no es fa de la cadena a guardar (apply chain) contra les coordenades de les cadenes ja guardades ->
#//de mateixa sequencia (mateix ordinal), sino contra TOTES les cadenes guardades
#//Es mes lent, pero he tingut problemes en aquest sentit amb l'ATP syntasa

#//2. Limitant chains: ara es poden fixar moltes limitant chains (no nomes una), cada una de les quals te el seu propi limit marcat a la opcio ->
#// --max_chains. Ara, ambdues opcions guarden dades en format llista, i aquestes llistes son combinades en el diccionari de dobles-llistes ->
#//Limitant_ordinals

#//3. Proportions multiplier option: Un multiplicador que multiplica els lles proporcions establertes a max_chains. Per exemple, si tenim una ->
#//estructura AABAABAABAAB que volem repetir 25 cops fiquem: --max_chains 1 2 --limitant_chains B A --proprtions_mulitplier 25

#//4. Seed option: Una altre opcio per donar als usuaris la possiblitat de triar el seed que vulguin 

#IMPORTANT: Soc conscient que només s'ensambla la meitat de l'atp syntasa. Ens falta una interacció que l'script de l'Alvaro no ha tret. 