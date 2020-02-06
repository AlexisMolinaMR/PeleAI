from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import ForwardSDMolSupplier
import pandas as pd
import warnings
import glob
import os
import time
import argparse as ap

start_time = time.time()

non_applying_descriptors = []

calc = Calculator(descriptors)

# path = "/home/alexis/Desktop/PeleAI_data/simulation_pdb/ligand_extract/ligands_2D/"

def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to extracted ligands.")
    
    args = parser.parse_args()

    path = args.input

    return path


def descriptor_calc(path):

    for filename in glob.glob(os.path.join(path, '*.sdf')):

        suppl = ForwardSDMolSupplier(filename)

        print("Calculating decriptors for " + os.path.basename(filename))
        df = calc.pandas(suppl)
        df.insert(0, 'Ligand', os.path.basename(filename).split(".")[0])

        if os.path.exists('dummy.xlsx'):
            append_write = 'a'
        else:
            append_write = 'w'

        with pd.ExcelWriter('dummy.xlsx', mode = append_write, engine="openpyxl") as writer:
            print("Writing descriptors")
            df.to_excel(writer)

    writer.save()

def descriptor_writer():

    descriptor_df = pd.read_excel('dummy.xlsx', None)
    all_ligand_descriptors = []

    for key in descriptor_df.keys():
        all_ligand_descriptors.append(descriptor_df[key])

    descriptors_concatenated = pd.concat(all_ligand_descriptors,axis=0,ignore_index=True)
    descriptors_concatenated = descriptors_concatenated.drop(columns=['Unnamed: 0'])

    for i in range(1,len(descriptors_concatenated.columns)):
        if (str(descriptors_concatenated.dtypes[i]) != 'float64') and (str(descriptors_concatenated.dtypes[i]) != 'int64') and (str(descriptors_concatenated.dtypes[i]) != 'bool'):
            non_applying_descriptors.append(descriptors_concatenated.columns[i])

    for descriptor in non_applying_descriptors:
        descriptors_concatenated = descriptors_concatenated.drop(columns=[descriptor])
        print("Dropping " + descriptor)

    writer = pd.ExcelWriter('final_descriptors_try.xlsx')
    descriptors_concatenated.to_excel(writer, sheet_name='descriptors', index=False)
    descriptors_concatenated.to_csv('final_descriptors_try.csv')
    writer.save()



def main():

        path = parseArg()

        descriptor_calc(path=path)

        descriptor_writer()

        print("DONE!")


if __name__ == "__main__":
    main()

    print("Descriptors calculated and stored in {} seconds.".format(time.time()-start_time))
