import sys
import os
import glob
import time

'''
import subprocess

bashCommand_start = "conda activate my-rdkit-env"
process = subprocess.Popen(bashCommand_start.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand_end = "conda deactivate"
process = subprocess.Popen(bashCommand_end.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

path = sys.argv[1] # path where PDB's are stored
save_path = sys.argv[2] # save png's with 2d ligand structure
ligand_pdb = sys.argv[3] # ligand name as in pdb file ## extend to set of ligands
keep_ligand = sys.argv[4] # yes/no ## remove or save ligands

for filename in glob.glob(os.path.join(path, '*.pdb')):
    ligand_name = (os.path.basename(filename).split("_trajectory")[0])
    file = open(filename).readlines()
    with open(save_path + ligand_name + ".pdb", "w") as ligand_file:
        for lines in file:
            k = lines.split()
            try:
                if k[3] == ligand_pdb:
                    ligand_file.write(lines)
            except:
                pass

    ligand = AllChem.MolFromPDBFile(save_path + ligand_name + ".pdb")
    Draw.MolToFile(ligand, save_path + ligand_name + '.png')

    if keep_ligand == "yes":
        continue
    else:
        os.remove(save_path + ligand_name + ".pdb")
