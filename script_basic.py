#!/usr/bin/python3

from Bio.PDB import *
from sys import *
import os

io = PDBIO()
os.system("rm temp*.pdb")

def savepdb (pdb_ob,filename):
	"""
	Save an object in pdb format
	"""
	io.set_structure(pdb_ob)
	io.save(filename)

def clashtest(chain1, chain2):
	"""
	Returns True if any of the atoms on chain1 clash with any of the atoms in chain2
	"""
	clashcounter = 0

	#Obtain a list of atoms for chain1 and a list of coordinates for chain 2
	atom1 = [ atom for atom in chain1.get_atoms() ]
	cordsatom2 = [ atom.get_coord() for atom in chain2.get_atoms() ]

	#Obtain number of atoms in both chains
	numatoms1 = len(atom1)
	numatoms2 = len(cordsatom2)

	#Prepare list of atoms for neighbour search
	ns = NeighborSearch(atom1)

	#Check collisions of every atom in chain2 to already-prepared-list of atoms1 (collision distance: 1 amstrong)
	for clashes in map(lambda x: ns.search(x,2.0), cordsatom2):
		#Count number of atoms clashing
		if len(clashes)> 0:
			clashcounter += 1

	#If any atom clashes, calculate the overall percentage of clashing atoms in chain1
	#//El missatge d'error el deixo de moment com a demostracio. A la versio definitiva no hi hauria de ser
	clashpercent  = clashcounter/numatoms1 * 100
	if clashcounter > 0:
		stderr.write("%i per cent atoms in chain %s are clashing with chain %s\n" %(clashpercent, chain1.get_id(), chain2.get_id()))

	#Return percentage of clashing atoms
	return clashpercent

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

	sup.apply(apply_chain)
	return apply_chain

def recursive(chain_ob, camefrom, interactions, allpdb):
	"""
	NOTA A MARTA: Aquesta subrutina, benvolguda Marta, potser ens solucioni gran part del treball o potser ens peti l'ordinador. L'important es que sortirem de dubtes
	L'he comentat molt perque la puguis entendre facil i rapidament
	ALERTA!!!!: si visualitzes els resultats a chimera es veuran malament la majoria dels cops, doncs aquest algoritme te tendencia a repetir subunitats. I si dues subunitats estan molt a prop, chimera les separa
	"""
	chainame = chain_ob.get_id()
	#//Tancament de seguretat per evitar que se li vagi l'olla
	if counter == 5:
		exit()
	#Iterate throught the chains which interact with the present chain
	for chainteracting in interactions[chainame]:
		#If this loop corresponds to the already superimposed chain where the present chain_ob comes from, skip it
		if chainteracting == camefrom:
			continue

		#print("chain fixed: %s\nchain to move: %s" % (chain_ob,chainteracting))
		#Extract the pdb where chain_ob and chainteracting are both toghether
		pdbname = interactions[chainame][chainteracting]
		pdb_ob = allpdb[pdbname]
		#chaintomove: chain of the same type than chain_ob in the pdb-file
		chaintomove = pdb_ob[0][chainame]
		#aplichain: chain that will be moved to interact with chain_ob
		applychain = pdb_ob[0][chainteracting]

		#Check collistions
		clashtest(chain_ob, applychain)

		#Superimpose
		appliedchain = superimpose(fix_chain = chain_ob, mov_chain = chaintomove, apply_chain = applychain)

		#Check if appliedchain collides with any of the 
		for chainmoved in chainsmoved[appliedchain.get_id()]:
			clashperc = clashtest(appliedchain, chainmoved)
			if clashperc > 80:
				return 
		#print("chain %s moved succesfully\n" % (appliedchain))
		#Add applied chain to coresponding list of already-moved chains
		chainsmoved[appliedchain.get_id()].append(appliedchain)
		#Add one to the superimpositions counter
		global counter
		counter += 1

		#Save appliedchain (the moved one) in a new pdb
		tempname = "temp" + str(counter) + ".pdb"
		savepdb(appliedchain, tempname)

		#Repeat process for appliedchain
		recursive(appliedchain, chainame, interactions, allpdb)


##########################################
##Store Relations between files and chains
##########################################

"""
Dictionary, list and other variables index:
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
	5. counter: counts the number of chains moved for bulding the present model
	6. Chainsmoved: dictionary with
		keys: all chains in the present model
		values: chain_objects moved of the corresponding chain type
"""

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

#declare some variables

#Create allpdb dictionary
allpdb = { filename : parser.get_structure(filename, filename) for filename in argv if argv.index(filename) != 0}

#Create subunits set and chainsbyfile dictionary
subunits = set()
chainsbyfile = {}
for pdb in allpdb:
	chainsbyfile[pdb] = list()
	for subunit in allpdb[pdb].get_chains():
		chainsbyfile[pdb].append(subunit)
		subunits.add(subunit.get_id())

#Create intersections dictionary
interactions = {}
intset = dict()

for element in subunits: 
	for filename in allpdb:
		chainames = [ a.get_id() for a in chainsbyfile[filename] ]
		if element in chainames:
			index = chainames.index(element)
			if index == 0:
				intset[chainames[1]] = filename
				
			else:
				intset[chainames[0]] = filename
	interactions[element] = intset
	intset = dict()

#Counter of chains
global counter
counter = 0

#Create chainsmoved dictionary
chainsmoved = { chain: [] for chain in subunits }

#print("allpdb: ",allpdb)
#print("subunits: ",subunits)
#print("chainsbyfile: ",chainsbyfile)
#print("interactions: ",interactions)
#print("chainsmoved: ",chainsmoved)

"""
//NOTA A SANDVITX: He canviat els pdbs que utiltzem com a probes pels que fan servir Oriol i Alvaro. Perque? perque no hi havia caigut en que, malgrat estiguin en fitxers separats, les nostres
cadenes ja estan positionades on els toca, d'ons venen d'un fitxer comu
Ara surten un munt de warnings, pero no t'alarmis, funciona igualment. Es perque els pdbs que utilitzen Alvaro i Oriol els falta una columna 
"""


"""
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
			pdbinteract = interactwith[pdb2]
			done = predone.copy()
			done.add(pdbinteract)
			if len(predone) == len(done): 
				continue
			else:
				count += 1
				testwhere = chainsbyfile[pdbinteract][option[0]].get_id()
				print(testwhere)
				if testwhere == fixchain:
					posfix = option[0]
					posmov = option[1]
				else:
					posfix = option[1]
					posmov = option[0]
				#print(pdbname, option[0], pdbinteract, posfix, pdbinteract, posmov)
				moved_chain = superimpose(chainsbyfile[pdbname][0], chainsbyfile[pdbinteract][posfix], chainsbyfile[pdbinteract][posmov])
				io.set_structure(moved_chain)
				tempname = "temp" + str(count) + ".pdb"
				io.save(tempname)
		predone = done
	os.system("cat temp*.pdb > krosis.pdb")
"""



#Obtain a random pdb interacting file
seed = list(allpdb.keys())[0]
#Iterate over pdb seed file chains
areclashes = False
for chain in chainsbyfile[seed]:

	#Check collisions for seed chains
	for chainmoved in chainsmoved[chain.get_id()]:
		clashperc = clashtest(chain, chainmoved)
		if clashperc > 80:
			areclashes = True 

	#save seed chains 
	if not areclashes:
		savepdb(chain, "temp" + str(counter) + ".pdb")
		chainsmoved[chain.get_id()].append(chain)

		recursive(chain, chainsbyfile[seed][1].get_id(), interactions, allpdb)

os.system("cat temp*.pdb > krosis.pdb")
