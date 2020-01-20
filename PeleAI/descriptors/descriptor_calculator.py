from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import ForwardSDMolSupplier
import pandas as pd
import warnings
import glob
import os

calc = Calculator(descriptors)

#print(len(calc.descriptors))

#mols = [Chemfor filename in glob.glob(os.path.join(path, '*.sdf')):

path = "/home/alexis/Desktop/PeleAI_data/simulation_pdb/ligand_extract/ligands_2D/"

descriptors_df = pd.DataFrame()

for filename in glob.glob(os.path.join(path, '*.sdf')):
    suppl = ForwardSDMolSupplier(filename)
    df = calc.pandas(suppl)
    descriptors_df.append(df)
#    for mol in suppl:
#        print(list(calc(mol))[0])
