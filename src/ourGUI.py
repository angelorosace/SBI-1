#!/usr/bin/python3

from tkinter import filedialog, messagebox
from tkinter import *

class Window(Frame):

	def __init__(self, master=None):
		Frame.__init__(self, master)                 
		self.master = master
		self.init_window()

	#Creation of init_window
	def init_window(self):

		# changing the title of our master widget      
		self.master.title("MatchProt GUI")

		# allowing the widget to take the full space of the root window
		#self.pack(fill=BOTH, expand=1)

		# creating a menu instance
		menu = Menu(self.master)
		self.master.config(menu=menu)    

		# create file object, add it on the menu and add command actions
		file = Menu(menu)
		file.add_command(label="Open", command=self.open_files)
		file.add_command(label="Output file", command=self.output_file)
		file.add_command(label="Output path", command=self.output_path)
		file.add_command(label="Temp", command=self.temp_file)
		file.add_command(label="Seed", command=self.seed_file)
		file.add_separator()
		file.add_command(label="Exit", command=self.client_exit)

		#added "file" to our menu
		menu.add_cascade(label="File", menu=file)

		# create file object, add it on the menu and add command actions
		options = Menu(menu)
		options.add_command(label="Verbose", command=self.verbose_option)
		options.add_command(label="Max chains", command=self.max_chains)
		options.add_command(label="Limitant chains", command=self.limitant_chains)
		options.add_command(label="Proportions multipliers", command=self.proportions_multipliers)
		options.add_command(label="Unique ids", command=self.unique_ids)

		#added "file" to our menu
		menu.add_cascade(label="Options", menu=options)

		# create file object, add it on the menu and add command actions
		helpme = Menu(menu)
		helpme.add_command(label="Info", command=self.help_info)
		helpme.add_command(label="Infiles", command=self.help_open)
		helpme.add_command(label="Outfiles", command=self.help_output)
		helpme.add_command(label="Outfile path", command=self.help_path)
		helpme.add_command(label="Verbose", command=self.help_verbose)
		helpme.add_command(label="Temp", command=self.help_temp)
		helpme.add_command(label="Seed", command=self.help_seed)
		helpme.add_command(label="Max chains", command=self.help_max)
		helpme.add_command(label="Limitant chains", command=self.help_limitant)
		helpme.add_command(label="Proportions chains", command=self.help_proportions)
		helpme.add_command(label="Unique ids", command=self.help_unique)

		#added "file" to our menu
		menu.add_cascade(label="Help", menu=helpme)

	def open_files(self):
		infiles = filedialog.askopenfilenames(title = "Select all interactions pdb files", filetypes = [("pdb files","*.pdb")]) 

	def output_file(self):
		Label(root, text="Output File Name:").grid(row=0)
		entry = Entry(root)
		entry.grid(row=0, column=1) 

		# Calling on_change when you press the return key
		entry.bind("<Return>", self.on_change_file)

	def on_change_file(self, e):
		outfile = e.widget.get()
		messagebox.showinfo("Output File", outfile)

	def output_path(self):
		Label(root, text="Output Path:").grid(row=2)
		entry = Entry(root)
		entry.grid(row=2, column=1) 

		# Calling on_change when you press the return key
		entry.bind("<Return>", self.on_change_path)

	def on_change_path(self, e):
		path = e.widget.get()
		messagebox.showinfo("Output Path", path)

	def temp_file(self):
		Label(root, text="Temp file:").grid(row=4)
		temp_save = BooleanVar()
		save = Radiobutton(root, text="Save temporary files", variable=temp_save, value=True, height=2, width=20)
		delate = Radiobutton(root, text="Delate temporary files", variable=temp_save, value=False, height=2, width=20)
		save.grid(row=4, column=1)
		delate.grid(row=5, column=1)
		go = Button(root, text="Save", command=lambda var=temp_save: self.sel_temp(var))
		go.grid(row=6, column=1)

	def sel_temp(self, var):
		temp = var.get()
		messagebox.showinfo("Temp", "Will be the temporary files saved? " + str(temp))

	def seed_file(self):
		Label(root, text="Seed file:").grid(row=7)
		entry = Entry(root)
		entry.grid(row=7, column=1) 

		# Calling on_change when you press the return key
		entry.bind("<Return>", self.on_change_seed)

	def on_change_seed(self, e):
		seed_pdb = e.widget.get()
		messagebox.showinfo("Seed file", seed_pdb)

	def client_exit(self):
		exit()

	def verbose_option(self):
		Label(root, text="Verbose:").grid(row=8)
		verbose = BooleanVar()
		verbose_yes = Radiobutton(root, text="Display verbose", variable=verbose, value=True, height=2, width=20)
		verbose_no = Radiobutton(root, text="Don't show log messages", variable=verbose, value=False, height=2, width=20)
		verbose_yes.grid(row=8, column=1)
		verbose_no.grid(row=9, column=1)
		go = Button(root, text="Save", command=lambda var=verbose: self.sel_verbose(var))
		go.grid(row=10, column=1)

	def sel_verbose(self, var):
		verbose = var.get()
		messagebox.showinfo("Verbose", "Is verbose activated? " + str(verbose))

	def max_chains(self):
		Label(root, text="Max chains: ").grid(row=11)
		chains = Entry(root)
		chains.grid(row=11, column=1) 

		# Calling on_change when you press the return key
		chains.bind("<Return>", self.on_change_chains)

	def on_change_chains(self, e):
		max_chains = e.widget.get()
		messagebox.showinfo("Max chains", max_chains)

	def limitant_chains(self):
		Label(root, text="Limitant chains: ").grid(row=12)
		limitant = Entry(root)
		limitant.grid(row=12, column=1)

		# Calling on_change when you press the return key
		limitant.bind("<Return>", self.on_change_limitant)

	def on_change_limitant(self, e):
		limitant_chains = e.widget.get()
		messagebox.showinfo("Limitant chains", limitant_chains)
		
	def proportions_multipliers(self):
		Label(root, text="Proportions multipliers: ").grid(row=13)
		pmult = Entry(root)
		pmult.grid(row=13, column=1)

		# Calling on_change when you press the return key
		pmult.bind("<Return>", self.on_change_pmult)

	def on_change_pmult(self, e):
		proportions_multipliers = e.widget.get()
		messagebox.showinfo("Proportions multipliers", proportions_multipliers)

	def unique_ids(self):
		Label(root, text="Unique ids: ").grid(row=14)
		unique_option = BooleanVar()
		unique_yes = Radiobutton(root, text="Non-unique: different", variable=unique_option, value=True, height=2, width=20)
		unique_no = Radiobutton(root, text="Non-unique: equal", variable=unique_option, value=False, height=2, width=20)
		unique_yes.grid(row=14, column=1)
		unique_no.grid(row=15, column=1)

		go = Button(root, text="Save", command=lambda var=unique_option: self.sel_unique(var))
		go.grid(row=16, column=1)

	def sel_unique(self, var):
		unique = var.get()
		messagebox.showinfo("Unique", "Treate non-unique sequence chains as different chains. " + str(unique))


	def help_info(self):
		messagebox.showinfo("Info", "Matchprot reconstructs protein complexes from its individual protein interactions")

	def help_open(self):
		messagebox.showinfo("Open", "<Mandatory> Input PDB interaction files. Every file has to contain a unique, one-to-one, protein interaction")
	
	def help_output(self):
		messagebox.showinfo("Output", "Output file name")

	def help_path(self):
		messagebox.showinfo("Output path", " Output directory")

	def help_verbose(self):
		messagebox.showinfo("Verbose", "Print log comments")

	def help_temp(self):
		messagebox.showinfo("Temp", "Don't delete temporary files. Each one of them contain a subunit of the complex")

	def help_seed(self):
		messagebox.showinfo("Seed", "Name of the pdb file from which the reconstruction will start. Default is random")

	def help_max(self):
		messagebox.showinfo("Max chains", "Maximum number of subunits for the final model. If limitant_chains is activated, list of absolute proportions for subunit in complex")

	def help_limitant(self):
		messagebox.showinfo("Limitant chins", "List of limitant chain identifiers for the  model. Each sequence specified in this list will be in the model equal or less times than the corresponding number in --limitant_chains list multiplied by the --proportions_multiplier option. If a chain has more than one name in your model, put just one of its identifiers.")

	def help_proportions(self):
		messagebox.showinfo("Proportions multipliers", "If --limitant_chains is activated, number to multiply the proportions by. Default is 1")

	def help_unique(self):
		messagebox.showinfo("Unique ids", "Set this option if your model has chains with non-unique sequences BUT you still want to treat them as different chains")


def main(): 
	global root 
	root = Tk()
	app = Window(root)
	root.mainloop()

if __name__ == '__main__':
	main()





'''
	def seed_file(self):
		print(self.infiles)
		if not self.infiles:
			messagebox.showarning("Warning", "You have to select interaction files before selecting the seed file!")
		else:
			scrollbar = Scrollbar(root)
			scrollbar.pack( side = RIGHT, fill = Y )
			mylist = Listbox(root, yscrollcommand = scrollbar.set)
			for pdb in self.infiles:
				mylist.insert(END, str(pdb))

			mylist.pack( side = LEFT, fill = BOTH )
			scrollbar.config( command = mylist.yview )
'''