# Creation of an objective function for the RNA folding problem

============
## Context 

For a given ribonucleotide chain , the RNA folding problem consist in finding the native fold among the large number of possible conformation.The native fold being the one with the lowest Gibbs free energy,the objective function should be an estimator of this energy

The script will   

* Train the objective function using interactomic distance distribution computed from a dataset of experimentally determined 3D structures

* Plot the estimated Gibbs free energy as a function of the interatomic distance

* Use the objective function to evaluate predicted structures from the RNA
puzzles dataset

 ---
 ###  Training script : Compute interactomic distance from a given set of PDB files : 
                     
* only C3 atoms are taken into account
* only  intrachain distances are considered
* only consider residues separed by 3 position on the sequence
                      
* compute observed frequencies ( range of 0 to 20 Ångström of observing two residues i and j separared by a distance bin r is calculated as follows :

$$ f_{i,j} ^{OBS}(r) = { N_{i,j}(r) \over N_{i,j} } $$
                      
* compute reference frequency , different residue types are  indistinct X                    
 
 $$ f_{X,X} ^{REF}(r) = { N_{X,X}(r) \over N_{X,X}} $$

* compute log-ratio of the two frequency (score of pseudo energy)

$$ u_{i,j}(r) = { -log \left( f _{i,j} ^{OBS}(r) \over f_{i,j} ^{REF}(r) \right) } $$




---
# Installation

Clone the repository :
```bash
$ git clone https://github.com/hamidsta/RNA-folding-problem.git
```
Install the required dependency

```bash
$ pip install -r dependency.txt
```

---
# Running the script

In order to run the code , make sure that utils.py and main.py are located in the SAME FOLDER.
Main.py allow an interactive usage for the utilisator. I advise to check --help to understand how to make right use of the script .
Script only take RNA STRUCTURE 

Run the script on the terminal ( or you can just double click on the file )

```bash
$ python main.py
```

---
# Getting Data
The data can be downloaded  from https://www.rcsb.org/ 
Make sure the data are stored in an empty folder 




