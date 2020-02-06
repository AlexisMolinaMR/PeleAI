import sys
import os
import glob
import time
import argparse as ap

'''
import subprocess

bashCommand_start = "conda activate my-rdkit-env"
process = subprocess.Popen(bashCommand_start.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand_end = "conda deactivate"
process = subprocess.Popen(bashCommand_end.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

path = sys.argv[1] # path where PDB's are stored
save_path = sys.argv[2] # save png's with 2d ligand structure
ligand_pdb = sys.argv[3] # ligand name as in pdb file ## extend to set of ligands
keep_ligand = sys.argv[4] # yes/no ## remove or save ligands
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to PDB containing protein-ligand complex.")
    parser.add_argument("-o", "--output", required=True, type=str, help="path to save ligand-only structures.")
    parser.add_argument("-l", "--ligand", required=True, type=str, nargs="*", help="ligand(s) name in the given pdb files.")
    parser.add_argument("-r", "--remove", required=False, type=str2bool, help="a boolean. Decide whether to keep or delete the ligand-only pdb.")

    args = parser.parse_args()

    path = args.input
    save_path = args.output
    ligand_pdb = args.ligand
    keep_ligand = args.remove

    return path, save_path, ligand_pdb, keep_ligand


def extract_ligand(path, save_path, ligand_pdb):
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

        if keep_ligand:
            continue
        else:
            os.remove(save_path + ligand_name + ".pdb")

def main():

        path, save_path, ligand_pdb, keep_ligand = parseArg()

        extract_ligand(path=path, save_path=save_path, ligand_pdb=ligand_pdb, keep_ligand=keep_ligand)

        print("2D ligand structures available here: " + save_path)


if __name__ == "__main__":
    main()
