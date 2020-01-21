from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import ForwardSDMolSupplier
import pandas as pd
import warnings
import glob
import os
import time

start_time = time.time()

calc = Calculator(descriptors)

path = "/home/alexis/Desktop/PeleAI_data/simulation_pdb/ligand_extract/ligands_2D/"

for filename in glob.glob(os.path.join(path, '*.sdf')):

    suppl = ForwardSDMolSupplier(filename)

    print("Calculating decriptors for " + os.path.basename(filename))
    df = calc.pandas(suppl)

    if os.path.exists('dummy.xlsx'):
        append_write = 'a'
    else:
        append_write = 'w'

    with pd.ExcelWriter('dummy.xlsx', mode = append_write, engine="openpyxl") as writer:
        print("Writing descriptors")
        df.to_excel(writer)

writer.save()

descriptor_df = pd.read_excel('dummy.xlsx', None)
all_ligand_descriptors = []
for key in descriptor_df.keys():
    all_ligand_descriptors.append(descriptor_df[key])
descriptors_concatenated = pd.concat(all_ligand_descriptors,axis=0,ignore_index=True)
writer = pd.ExcelWriter('final_descriptors.xlsx')
descriptors_concatenated.to_excel(writer, sheet_name='descriptors', index=False)
writer.save()

print("Descriptors calculated and stored in {} seconds.".format(time.time()-start_time))
