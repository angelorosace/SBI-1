# Protein assembler

Proteins are not isolated on the cell. They interact with each other.
Some of this interactions are obligated, and result in two or more peptide chains forming an oligomer or **protein complex**.

The objective of protein assembler is to use individual, obligated interactions between different subunits to assemble the whole protein complex. 

## Installation

First download a zipped clone from git repository, either by browser or command line:

```
wget https://github.com/comradeameba/SBI/archive/master.zip
```

Uncompress the zip clone:

```
unzip master.zip
```

Once uncompressed, from you current directory:

1. **From binary:** Just add SBI-master/bin directory to your path. Alternatively, you can copy SBI-master/bin/protein_assembler to your path folder.
2. **From source script:** Just use the script SBI-master/src/protein_assembler.py. It's a little bit faster than the binary, but it has some extra requirements
  + python interpreter (version 3.0 or higher)
  + Biopython python module

## Usage

###Examples

For begginers:
```
protein_assembler -i pdb/hemoglobin/*
```

Intermediate user:
```
protein_assembler -i pdb/proteosome/* -o proteosome_example -d ./results/ -v
```

Advanced user:
```
protein_assembler -i pdb/2f1d_fosfate_dehydratase/* -v -o fosfate_dehydratase_example -d ./results/ -s pdb/2f1d_fosfate_dehydratase/XG.pdb -t -m 24  
```

Professional:
```
protein_assembler -i pdb/5ara_atp_syntasa/* -v -d ./results/ -o atp_syntase_example  -m 1 1 8 1 -l H W P T -s pdb/5ara_atp_syntasa/PH.pdb
```

### Options:

**-h/--help:** show help message and exit.

**-i/--input:** list of input pdb files. For this program to work properly, two different chains CAN NOT have the same id, even if they are in separate files. 

**-o/--output:** name of the output pdb file.

**-d/--path:** path where to store output pdb and temporary files

**-v/--verbose**: show log messages on command line.

**-s/--seed**: select an input seed pdb. The seed pdb is the one from which the protein complex assembling will start. Default is random.

**-t/--temp**: conserve temporary files. Protein_assembler will create a pdb file for every subunit saved for the final model, but once the execution all temporary files will be removed. Use this option if you want to preserve this temporary pdb files.

**-m/--max_chains:**  it has two forms of usage: 

1. As a single integer, --max_chains is the maximum number allowed of subunits for this protein complex. The program will stop when the limit is reached.
2. As a list of integers, --max_chains are the limits for the subunits specified in -l option (see atp syntase example on analysis_examples file  for more information about this option).

**-l/--limitant_chains:** if -m is a list of integers, -l should contain a list of same length with ids from desired chains. This way, the maximum number of each one of this chains is limited by the corresponding integer in -m option. If a chain has different ids on different pdbs, just put one of them (the program will apply the limit to all chains with same sequence). Usefull for assymetrical structues (see atp syntase example on analysis_examples file for more information about this option).

**-u/--unique_ids:** If this option is activated, all input chains with different chain ids will be treated as different subunits, regardless of the sequence composition (virus example on analysis_examples file for more information about this option).

**-f/--limitant_chains_multiplier:** Multiply all the max_chains integers by an integer. Usefull for large, repetitive structures, like microtubules.
  
## Citation

  + **Biopython module:** Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423

