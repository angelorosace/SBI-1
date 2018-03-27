#!/usr/bin/python3

from Bio.PDB import *
from sys import *
import os
from tkinter import filedialog, constants
from tkinter import *

root = Tk()
filenames =  filedialog.askopenfilenames(title = "Select all interactions pdb files", filetypes = [("pdb files","*.pdb")]) 

# Create a PDBparser object
parser = PDBParser(PERMISSIVE = 1)

#Create allpdb dictionary. THIS HAS TO BE CHANGED IN THE ORIGINAL SCRIP IF THE PROGRAM IS NOT EXECUTED BY COMMAND LINE.
if len(argv) == 1: 
	allpdb = { filename : parser.get_structure(filename, filename) for filename in filenames }
else:
	print("We will use command line version")
	
#print(allpdb)
"""
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

#print(chainsbyfile)



toplevel = Tk()
toplevel.withdraw()

filename_list = list()
filename = filedialog.askopenfilename(title = "Select file",filetypes = (("pdb files","*.pdb"),("all files","*.*")))
if os.path.isfile(filename):
	filename_list.append(filename)
else: 
	print('No file chosen')
	raw_input('Ready, push Enter')

print(filename_list)


def open_file(self, whatever = None, filename = None):
        if not filename:
            self.filename = tkinter.filedialog.askopenfilename(filetypes = self._filetypes)
        else:
            self.filename = filename
        if not (self.filename == ''):
            f = open(self.filename, 'r')
            f2 = f.read()
            self.text.delete('1.0', 'end')
            self.text.insert('1.0', f2)
            f.close()
            self.text.title('Sticky %s)' % self.filename)




import tkinter
import tkinter.messagebox

top = tkinter.Tk()

def showMessage():
        tkinter.messagebox.showinfo("Wellcome!", "Hello PYT students!")

b = tkinter.Button( top, text="Say Hello", command=showMessage)

b.pack()

top.mainloop()
"""