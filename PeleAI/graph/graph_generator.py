import os
import os.path
import glob
import math
import sys
import csv
import time
import argparse as ap
import numpy as np
import pandas as pd
from prody import *


def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--pose", required=True, type=str, help="path to docking poses")
    parser.add_argument("-l", "--ligand", required=True, type=str, help="ligand name")
    parser.add_argument("-r", "--radius", required=True, type=float, help="binding pocket selection radius")
    parser.add_argument("-at", "--atom_types", required=False, action='store_true', help="compute matrices for atom types instead of elements")
    parser.add_argument("-gc", "--geomcenter", required=False, action='store_true', help="if selected, method will proceed with geometric ligand center" )
    parser.add_argument("-ds", "--doubleselection", required=False, action='store_true', help="if selected, method will proceed with two selection spheres at outmost atoms of ligand's x-axis" )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-ex', "--exponential", action='store_true', help="use exponential as decay function")
    group.add_argument('-lo', "--lorentz", action='store_true', help="use Lorentz functions as decay function")
    parser.add_argument('-lp', "--inligprot", required=False, action='store_true', help="include intraligand and intraprotein interactions")
    parser.add_argument('-n', '--file_name', required=False, type=str, help='custom name to output csv files')

    args = parser.parse_args()

    path = args.pose
    ligand_name = args.ligand
    selection_radius = args.radius
    at_types = args.atom_types
    if args.geomcenter is not None:
        geom_center = args.geomcenter
    else:
        geom_center = False
    if args.doubleselection is not None:
        double_center = args.doubleselection
    else:
        double_center = False
    if args.exponential is not None:
        expo = args.exponential
        lor = None
    else:
        lor = args.lorentz
        expo = None
    if args.inligprot is not None:
        inligprot = args.inligprot
    else:
        inligprot = None

    if args.file_name is not None:
        name = args.file_name + '_'
    else:
        name = ""


    return path, ligand_name, at_types, selection_radius, geom_center, double_center, expo, lor, inligprot, name



def PDBParser(path):
    '''
    This function will parse the different docking poses
    in PDB format.
    '''

    pose_store = {}

    with open(path, 'r') as pdb:
        for lines in pdb:
            line = lines.strip()
            if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20].strip() in ["MSE","KCX"]):
                if not pose_store:
                    pose_store[path.split("/")[-1].split(".")[0]] = line + '\n'
                else:
                    pose_store[path.split("/")[-1].split(".")[0]] += line + '\n'
            elif line[0:6] == "HETATM":
                if not pose_store:
                    pose_store[path.split("/")[-1].split(".")[0]] = line + '\n'
                else:
                    pose_store[path.split("/")[-1].split(".")[0]] += line + '\n'

    p = parsePDB(path)

    return pose_store, p

def binding_pocket_selection(pose_store, p, ligand_name, selection_radius, center, double_center):
    '''
    This function will find by default mass center
    of the ligand using Prody.
    If the -gc option is selected the spatial center of
    the ligand is used by computing the mean distance between the
    furthest x axis coordinates.
    If the -ds option is selected a dobule sphere procedure is followed
    by selecting the furthest x axis atoms.
    '''

    amino = ['CYS', 'ASP', 'SER', 'GLN', 'LYS','ILE', 'PRO', 'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
    two_let_atom_code = ['Br','FE']

    coord = []
    x_coord = []
    binding_pocket = []
    ligand = []

    min_coord = None
    max_coord = None


    for pose in pose_store:
        structure = pose_store[pose].split('\n')

    for line in structure:
        if line[17:20] == ligand_name:
            coord.append(float(line[30:38]))
            coord.append(float(line[38:46]))
            coord.append(float(line[46:54]))

            ligand_atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-2:].strip())
            ligand.append(ligand_atom)
            ligand_atom = ()


    for i in range(0, len(coord), 3):
        x_coord.append(float(coord[i]))



    if double_center:
        print('\nDouble center of ligand selected')

        x_out_left = coord[coord.index(min(x_coord))]
        y_out_left = coord[coord.index(min(x_coord))+1]
        z_out_left = coord[coord.index(min(x_coord))+2]
        x_out_right = coord[coord.index(max(x_coord))]
        y_out_right = coord[coord.index(max(x_coord))+1]
        z_out_right = coord[coord.index(max(x_coord))+2]

        print('\n')
        print("Left sphere center coordinates: ", x_out_left, y_out_left, z_out_left)
        print("Right sphere center coordinates: ", x_out_right, y_out_right, z_out_right)
        print("Spheres radii: ", selection_radius)

        for line in structure:
            if (line[0:6].strip() == "ATOM" or line[0:6].strip() == "HETATM") and line[17:20].strip() in amino:

                x1 = math.pow((float(line[30:38]) - x_out_left), 2)
                y1 = math.pow((float(line[38:46]) - y_out_left), 2)
                z1 = math.pow((float(line[46:54]) - z_out_left), 2)

                if (x1 + y1 + z1) <= selection_radius**2:
                    atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-2:].strip())
                    binding_pocket.append(atom)
                    atom = ()

        for line in structure:
            if (line[0:6].strip() == "ATOM" or line[0:6].strip() == "HETATM") and line[17:20].strip() in amino:

                x1 = math.pow((float(line[30:38]) - x_out_right), 2)
                y1 = math.pow((float(line[38:46]) - y_out_right), 2)
                z1 = math.pow((float(line[46:54]) - z_out_right), 2)

                if (x1 + y1 + z1) <= selection_radius**2:
                    if line[-3:].strip() in two_let_atom_code:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip())

                    else:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip()[0])
                    binding_pocket.append(atom)
                    atom = ()

        print("Number of protein atoms selected: {}".format(len(binding_pocket)))
        print('Number of ligand atoms selected: {}'.format(len(ligand)))
        print('Total number of atoms selected: {}'.format(len(binding_pocket)+len(ligand)))

    elif center:

        print("\nGeometric ligand center selected")

        sphere_center_x = (coord[coord.index(max(x_coord))] + coord[coord.index(min(x_coord))]) / 2
        sphere_center_y = (coord[coord.index(max(x_coord))+1] + coord[coord.index(min(x_coord))+1]) / 2
        sphere_center_z = (coord[coord.index(max(x_coord))+2] + coord[coord.index(min(x_coord))+2]) / 2

        print('\n')
        print("Sphere center coordinates: ", sphere_center_x, sphere_center_y, sphere_center_z)
        print("Sphere radius: ", selection_radius)

        for line in structure:
            if (line[0:6].strip() == "ATOM" or line[0:6].strip() == "HETATM") and line[17:20].strip() in amino:

                x1 = math.pow((float(line[30:38]) - sphere_center_x), 2)
                y1 = math.pow((float(line[38:46]) - sphere_center_y), 2)
                z1 = math.pow((float(line[46:54]) - sphere_center_z), 2)

                if (x1 + y1 + z1) <= selection_radius**2:
                    if line[-3:].strip() in two_let_atom_code:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip())

                    else:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip()[0])
                    binding_pocket.append(atom)
                    atom = ()

        print("Number of protein atoms selected: {}".format(len(binding_pocket)))
        print('Number of ligand atoms selected: {}'.format(len(ligand)))
        print('Total number of atoms selected: {}'.format(len(binding_pocket)+len(ligand)))

    else:
        print('\nLigand mass center selected')
        ligand_selection = p.select('not water and hetero')

        weights = ligand_selection.getMasses()
        mass_center = calcCenter(ligand_selection, weights)


        sphere_center_x = mass_center[0]
        sphere_center_y = mass_center[1]
        sphere_center_z = mass_center[2]

        print("\nSphere center coordinates: {}, {}, {}".format(sphere_center_x, sphere_center_y, sphere_center_z))
        print("Sphere radius: {}".format(selection_radius))

        for line in structure:
            if (line[0:6].strip() == "ATOM" or line[0:6].strip() == "HETATM") and line[17:20].strip() in amino:

                x1 = math.pow((float(line[30:38]) - sphere_center_x), 2)
                y1 = math.pow((float(line[38:46]) - sphere_center_y), 2)
                z1 = math.pow((float(line[46:54]) - sphere_center_z), 2)

                if (x1 + y1 + z1) <= selection_radius**2:
                    if line[-3:].strip() in two_let_atom_code:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip())

                    else:
                        atom = (line[17:20].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip(), line[-3:].strip()[0])
                    binding_pocket.append(atom)
                    atom = ()

        print("Number of protein atoms selected: {}".format(len(binding_pocket)))
        print('Number of ligand atoms selected: {}'.format(len(ligand)))
        print('Total number of atoms selected: {}'.format(len(binding_pocket)+len(ligand)))


    print(binding_pocket)

    return binding_pocket, ligand



def elementsDistanceCalc(binding_pocket, ligand, inligprot):
    '''
    Compute distances between selected elements and discard covalent interactions.
    '''

    elements = []
    at_count = 0

    atom_interactions = []

    used_elements = []

    atomic_radius_elements = {'C':0.7, 'O':0.6, 'N':0.65, 'H':0.25, 'Cl':1, 'P':1, 'S':1, 'FE':1.4, 'Br':1.15, 'F':0.5, 'I':1.4}

    if inligprot:
        for atom in binding_pocket:
            if atom not in elements:
                atom_types.append(atom[-1])
        for atom in ligand:
            if atom not in elements:
                elements.append(atom[-1])

        elements = list(set(elements))

        for element in elements:
            if (len(element) > 2) and (element[0] in atom_types):
                elements.remove(element)

        print("Number of possible subgraphs: {}".format(((len(elements)**2)+len(elements))/2))
        print('\n')

        protein_ligand = binding_pocket + ligand

        protein_ligand.sort(key=lambda tup: tup[-1])

        elements.sort()

        for i in range(len(protein_ligand)):
            try:
                if protein_ligand[i][-1] == elements[at_count]:

                    used_elements.append(elements[at_count])
                    at_count += 1

                    for j in range(i, len(protein_ligand)):
                        if i != j:

                            x1 = math.pow(float(protein_ligand[j][2]) - float(protein_ligand[i][2]), 2)
                            y1 = math.pow(float(protein_ligand[j][3]) - float(protein_ligand[i][3]), 2)
                            z1 = math.pow(float(protein_ligand[j][4]) - float(protein_ligand[i][4]), 2)

                            print(protein_ligand[i][-1])
                            if x1 + y1 + z1 <= atomic_radius_elements[protein_ligand[i][-1][-3:]]**2:
                                print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                            else:
                                distance_interaction = (protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1))
                                atom_interactions.append(distance_interaction)
                                distance_interaction = ()
            except:
                pass


        ligand_elements = None
        protein_elements = None


    else:

        at_count = 0

        protein_elements = ['C','O','N','S']
        protein_elements.sort()
        ligand_elements = ['H','C','N','O','F','P','S','Cl','Br','I']
        ligand_elements.sort()

        ligand.sort(key=lambda tup: tup[-1])
        binding_pocket.sort(key=lambda tup: tup[-1])

        for i in range(len(ligand)):
            if binding_pocket[i][-1] == protein_elements[at_count]:
                at_count += 1
            for j in range(i, len(ligand)):

                x1 = math.pow(float(ligand[j][2]) - float(binding_pocket[i][2]), 2)
                y1 = math.pow(float(ligand[j][3]) - float(binding_pocket[i][3]), 2)
                z1 = math.pow(float(ligand[j][4]) - float(binding_pocket[i][4]), 2)

                try:
                    if int(ligand[i][-1][0]):

                        if (x1 + y1 + z1 <= atomic_radius_elements[binding_pocket[i][-1]]**2) or (x1 + y1 + z1 <= atomic_radius_elements[ligand[i][1][0]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()
                except:
                    if (x1 + y1 + z1 <= atomic_radius_elements[binding_pocket[i][-1]]**2) or (x1 + y1 + z1 <= atomic_radius_elements[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()


        elements = None

    return atom_interactions, elements, ligand_elements, protein_elements


def atomTypesDistanceCalc(binding_pocket, ligand):

    atom_interactions = []

    atom_types = []

    used_types = []

    atomic_radius_types = {'C':0.7,'CA':0.7,'CB':0.7,'CG':0.7,'CD':0.7,'CE':0.7,'CZ':0.7,'CH':0.7, 'O':0.6, 'N':0.65, 'H':0.25, 'Cl':1, 'P':1, 'S':1, 'FE':1.4, 'Br':1.15, 'F':0.5, 'I':1.4}

    at_count = 0

    protein_atom_types = ['C','CA','CB','CG','CD','CE','CZ','CH','O','N','S']
    protein_atom_types.sort()
    ligand_atom_types = ['H','C','N','O','F','P','S','Cl','Br','I']
    ligand_atom_types.sort()

    ligand.sort(key=lambda tup: tup[-1])
    binding_pocket.sort(key=lambda tup: tup[-1])

    for i in range(len(ligand)):
        if binding_pocket[i][1] == protein_atom_types[at_count]:
            at_count += 1
        for j in range(i, len(ligand)):

            x1 = math.pow(float(ligand[j][2]) - float(binding_pocket[i][2]), 2)
            y1 = math.pow(float(ligand[j][3]) - float(binding_pocket[i][3]), 2)
            z1 = math.pow(float(ligand[j][4]) - float(binding_pocket[i][4]), 2)

            try:
                if int(ligand[i][-1][0]):

                    if (x1 + y1 + z1 <= atomic_radius_types[binding_pocket[i][-1]]**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][1][0]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()
            except:
                if (x1 + y1 + z1 <= atomic_radius_types[binding_pocket[i][-1]]**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                    print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1],protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                else:
                    distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                    atom_interactions.append(distance_interaction)
                    distance_interaction = ()

    for i in sorted(protein_atom_types):
        for j in sorted(ligand_atom_types):
            atom_types.append(j + '-' + i)

    atom_types.sort()


    return atom_interactions, atom_types, ligand_atom_types, protein_atom_types



def subgraphsWeights(atom_interactions, types, expo, lorentz, inligprot, ligand_atom_types, protein_atom_types, at_types):
    '''
    Compute subgraph weights using either exponential or Lorentz decay functions.
    '''

    weights = []
    count = 0
    atom_atom = []
    atom_rev = []
    uni_inter = []

    final_weigths = []

    if at_types:
        if expo is not None:
            for inter in atom_interactions:
                type_weight = (np.exp(-abs(inter[-1])/3), "{}-{}".format(inter[2],inter[1]))
                weights.append(type_weight)
                type_weight = ()

        elif lor is not None:
            for inter in atom_interactions:
                type_weight = ((1/1 + (abs(inter[-1])/3)**5), "{}-{}".format(inter[2],inter[1]))
                weights.append(type_weight)
                type_weight = ()

    else:
        if expo is not None:
            for inter in atom_interactions:
                type_weight = (np.exp(-abs(inter[-1])/3), "{}-{}".format(inter[0],inter[2]))
                weights.append(type_weight)
                type_weight = ()

        elif lor is not None:
            for inter in atom_interactions:
                type_weight = ((1/1 + (abs(inter[-1])/3)**5), "{}-{}".format(inter[0],inter[2]))
                weights.append(type_weight)
                type_weight = ()


    if inligprot:

        for i in range(len(types)):
            for j in range(len(types)):
                atom_atom.append(("{}-{}".format(types[i], types[j])))
                if types[i] != types[j] and "{}-{}".format(types[i], types[j]) not in atom_rev:
                    atom_rev.append(("{}-{}".format(types[j], types[i])))

        for ele in atom_rev:
            if ele in atom_atom:
                atom_atom.remove(ele)

        weights = list(set(weights))

        for w in weights:
            if w[-1] in atom_atom:
                final_weigths.append(w)

    else:

        for i in range(len(ligand_atom_types)):
            for j in range(len(protein_atom_types)):
                atom_atom.append(("{}-{}".format(ligand_atom_types[i], protein_atom_types[j])))
                if ligand_atom_types[i] != protein_atom_types[j] and "{}-{}".format(ligand_atom_types[i], protein_atom_types[j]) not in atom_rev:
                    atom_rev.append(("{}-{}".format(protein_atom_types[j], ligand_atom_types[i])))


        weights = list(set(weights))

        if not at_types:
            for w in weights:
                if w[-1] in atom_atom:
                    final_weigths.append(w)

        else:

            final_weigths = weights

    atom_atom.sort()

    return final_weigths, atom_atom


def laplacianMatrix(weights, atom_combinations, at_types):
    '''
    Compute Laplacian matrices.
    '''

    laplacian_matrices = []

    for element in atom_combinations:
        element_count = 0
        for weight in weights:
            if weight[-1] == element:
                element_count += 1
        print("\nLaplacian matrix for subgraph G({}). Matrix size {}x{}:\n".format(element, element_count, element_count))
        lp_mat = np.zeros(shape=(element_count, element_count))
        L_count = element_count - 1
        sum_count = 0

        for i in range(element_count):
            sum = 0
            for j in range(element_count):
                if j == L_count:
                    lp_mat[i][j] == 0
                elif j == sum_count:
                    for k in range(len(weights)):
                        if k != L_count:
                            sum += weights[k][0]
                    lp_mat[i][j] = sum
                else:
                    lp_mat[i][j] = -weights[j][0]

            L_count -= 1
            sum_count += 1

        print(lp_mat)

        subgraph_matrix = (element, lp_mat)
        laplacian_matrices.append(subgraph_matrix)
        subgraph_matrix = ()


    return laplacian_matrices

def adjacencyMatrix(weights, atom_combinations, at_types):
    '''
    Compute adjacency matrices.
    '''

    adjacency_matrices = []


    for element in atom_combinations:
        element_count = 0
        for weight in weights:
            if weight[-1] == element:
                element_count += 1
        print("\nAdjacency matrix for subgraph G({}). Matrix size {}x{}:\n".format(element, element_count, element_count))
        ad_mat = np.zeros(shape=(element_count, element_count))
        L_count = element_count - 1
        sum_count = 0

        for i in range(element_count):
            for j in range(element_count):
                if j == L_count:
                    ad_mat[i][j] == 0
                elif j == sum_count:
                    ad_mat[i][j] = 0
                else:
                    ad_mat[i][j] = -weights[j][0]

            L_count -= 1
            sum_count += 1

        print(ad_mat)

        subgraph_matrix = (element, ad_mat)
        adjacency_matrices.append(subgraph_matrix)
        subgraph_matrix = ()


    return adjacency_matrices


def laplacianStats(matrices,path):
    '''
    Compute statistics out of Laplacian matrices
    '''

    laplacian_statistics = []

    print("\nComputing statistics for laplacian matrices...")

    laplacian_statistics.append(path.split("/")[-1].split(".")[0])

    for mat in matrices:
        if mat[-1].size != 0:
            print('\n')
            print(mat[0])
            print('\n')
            print(np.linalg.eig(mat[-1])[1])
            maximum = np.linalg.eig(mat[-1])[1].max()
            print("\nMax eigenvalue {}".format(maximum))
            minimum = np.linalg.eig(mat[-1])[1].min()
            print("Min eigenvalue {}".format(minimum))
            sum_of = np.linalg.eig(mat[-1])[1].sum()
            print("Sum of eigenvalues {}".format(sum_of))
            mean = np.linalg.eig(mat[-1])[1].mean()
            print("Mean of eigenvalues {}".format(mean))
            median = np.median(np.linalg.eig(mat[-1])[1])
            print("Median of eigenvalues {}".format(median))
            std = np.linalg.eig(mat[-1])[1].std()
            print("Standard deviation of eigenvalues {}".format(std))
            vari = np.linalg.eig(mat[-1])[1].var()
            print("Variance of eigenvalues {}".format(vari))
            num = len(np.linalg.eig(mat[-1])[1])
            print("Number of eigenvalues {}".format(num))
            sec_pow = np.power(np.linalg.eig(mat[-1])[1], 2)
            total_sec_pow =  0
            for i in sec_pow:
                total_sec_pow += sum(i)
            print("Second power of eigenvalues {}".format(total_sec_pow))
        else:
            print('\n')
            print(mat[0])
            print('\n')
            maximum = 0
            print("\nMax eigenvalue {}".format(maximum))
            minimum = 0
            print("Min eigenvalue {}".format(minimum))
            sum_of = 0
            print("Sum of eigenvalues {}".format(sum_of))
            mean = 0
            print("Mean of eigenvalues {}".format(mean))
            median = 0
            print("Median of eigenvalues {}".format(median))
            std = 0
            print("Standard deviation of eigenvalues {}".format(std))
            vari = 0
            print("Variance of eigenvalues {}".format(vari))
            num = 0
            print("Number of eigenvalues {}".format(num))
            total_sec_pow = 0
            print("Second power of eigenvalues {}".format(total_sec_pow))

        laplacian_statistics.append(maximum)
        laplacian_statistics.append(minimum)
        laplacian_statistics.append(sum_of)
        laplacian_statistics.append(mean)
        laplacian_statistics.append(median)
        laplacian_statistics.append(std)
        laplacian_statistics.append(vari)
        laplacian_statistics.append(num)
        laplacian_statistics.append(total_sec_pow)

    return laplacian_statistics

def adjacencyStats(matrices, path):
    '''
    Compute statistics out of adjacency matrices.
    '''

    adjacency_statistics = []

    print("\nComputing statistics for adjacency matrices...")

    adjacency_statistics.append(path.split("/")[-1].split(".")[0])

    for mat in matrices:
        if mat[-1].size != 0:
            print('\n')
            print(mat[0])
            print('\n')
            print(np.linalg.eig(mat[-1])[1])
            maximum = np.linalg.eig(mat[-1])[1].max()
            print("\nMax eigenvalue {}".format(maximum))
            minimum = np.linalg.eig(mat[-1])[1].min()
            print("Min eigenvalue {}".format(minimum))
            sum_of = np.linalg.eig(mat[-1])[1].sum()
            print("Sum of eigenvalues {}".format(sum_of))
            mean = np.linalg.eig(mat[-1])[1].mean()
            print("Mean of eigenvalues {}".format(mean))
            median = np.median(np.linalg.eig(mat[-1])[1])
            print("Median of eigenvalues {}".format(median))
            std = np.linalg.eig(mat[-1])[1].std()
            print("Standard deviation of eigenvalues {}".format(std))
            vari = np.linalg.eig(mat[-1])[1].var()
            print("Variance of eigenvalues {}".format(vari))
            num = len(np.linalg.eig(mat[-1])[1])
            print("Number of eigenvalues {}".format(num))
            sec_pow = np.power(np.linalg.eig(mat[-1])[1], 2).sum()
            print("Second power of eigenvalues {}".format(sec_pow))

            if np.iscomplex(maximum):
                maximum = np.real(maximum, tol=1)
            adjacency_statistics.append(maximum)
            if np.iscomplex(minimum):
                minimum = np.real(minimum)
            adjacency_statistics.append(minimum)
            if np.iscomplex(sum_of):
                sum_of = np.real(sum_of)
            adjacency_statistics.append(sum_of)
            if np.iscomplex(mean):
                mean = np.real(mean)
            adjacency_statistics.append(mean)
            if np.iscomplex(median):
                median = np.real(median)
            adjacency_statistics.append(median)
            if np.iscomplex(std):
                std = np.real(std)
            adjacency_statistics.append(std)
            if np.iscomplex(vari):
                vari = np.real(vari)
            adjacency_statistics.append(vari)
            adjacency_statistics.append(num)
            adjacency_statistics.append(sec_pow)
        else:
            print('\n')
            print(mat[0])
            print('\n')
            maximum = 0
            print("\nMax eigenvalue {}".format(maximum))
            minimum = 0
            print("Min eigenvalue {}".format(minimum))
            sum_of = 0
            print("Sum of eigenvalues {}".format(sum_of))
            mean = 0
            print("Mean of eigenvalues {}".format(mean))
            median = 0
            print("Median of eigenvalues {}".format(median))
            std = 0
            print("Standard deviation of eigenvalues {}".format(std))
            vari = 0
            print("Variance of eigenvalues {}".format(vari))
            num = 0
            print("Number of eigenvalues {}".format(num))
            sec_pow = 0
            print("Second power of eigenvalues {}".format(sec_pow))

        adjacency_statistics.append(maximum)
        adjacency_statistics.append(minimum)
        adjacency_statistics.append(sum_of)
        adjacency_statistics.append(mean)
        adjacency_statistics.append(median)
        adjacency_statistics.append(std)
        adjacency_statistics.append(vari)
        adjacency_statistics.append(num)
        adjacency_statistics.append(sec_pow)

    return adjacency_statistics

def statsToCsv(laplacian_statistics, adjacency_statistics, expo, lorentz, inligprot, at_types, name):
    '''
    Write previously computed statistics to a csv file.
    '''

    if at_types:

        element_element = ['Br-CA', 'Br-CB', 'Br-CD', 'Br-CE', 'Br-CG', 'Br-CH', 'Br-CZ', 'Br-N', 'Br-O', 'Br-S', 'C-CA', 'C-CB', 'C-CD', 'C-CE', 'C-CG', 'C-CH', 'C-CZ', 'C-N', 'C-O', 'C-S', 'Cl-CA', 'Cl-CB', 'Cl-CD', 'Cl-CE', 'Cl-CG', 'Cl-CH', 'Cl-CZ', 'Cl-N', 'Cl-O', 'Cl-S', 'F-CA', 'F-CB', 'F-CD', 'F-CE', 'F-CG', 'F-CH', 'F-CZ', 'F-N', 'F-O', 'F-S', 'H-CA', 'H-CB', 'H-CD', 'H-CE', 'H-CG', 'H-CH', 'H-CZ', 'H-N', 'H-O', 'H-S', 'I-CA', 'I-CB', 'I-CD', 'I-CE', 'I-CG', 'I-CH', 'I-CZ', 'I-N', 'I-O', 'I-S', 'N-CA', 'N-CB', 'N-CD', 'N-CE', 'N-CG', 'N-CH', 'N-CZ', 'N-N', 'N-O', 'N-S', 'O-CA', 'O-CB', 'O-CD', 'O-CE', 'O-CG', 'O-CH', 'O-CZ', 'O-N', 'O-O', 'O-S', 'P-CA', 'P-CB', 'P-CD', 'P-CE', 'P-CG', 'P-CH', 'P-CZ', 'P-N', 'P-O', 'P-S', 'S-CA', 'S-CB', 'S-CD', 'S-CE', 'S-CG', 'S-CH', 'S-CZ', 'S-N', 'S-O', 'S-S']

        element_features = ['ligand']
        for elements in element_element:
            element_features.append(elements + "_max")
            element_features.append(elements + "_min")
            element_features.append(elements + "_sum")
            element_features.append(elements + "_mean")
            element_features.append(elements + "_median")
            element_features.append(elements + "_std")
            element_features.append(elements + "_vari")
            element_features.append(elements + "_num")
            element_features.append(elements + "_2dpow")

    else:

        element_element = ['Br-C', 'Br-N', 'Br-O', 'Br-S', 'C-C', 'C-N', 'C-O', 'C-S', 'Cl-C', 'Cl-N', 'Cl-O', 'Cl-S', 'F-C', 'F-N', 'F-O', 'F-S', 'H-C', 'H-N', 'H-O', 'H-S', 'I-C', 'I-N', 'I-O', 'I-S', 'N-C', 'N-N', 'N-O', 'N-S', 'O-C', 'O-N', 'O-O', 'O-S', 'P-C', 'P-N', 'P-O', 'P-S', 'S-C', 'S-N', 'S-O', 'S-S']
        element_features = ['ligand']
        for elements in element_element:
            element_features.append(elements + "_max")
            element_features.append(elements + "_min")
            element_features.append(elements + "_sum")
            element_features.append(elements + "_mean")
            element_features.append(elements + "_median")
            element_features.append(elements + "_std")
            element_features.append(elements + "_vari")
            element_features.append(elements + "_num")
            element_features.append(elements + "_2dpow")

    laplacian_st = []
    adjacency_st = []

    laplacian_st.append(laplacian_statistics)
    adjacency_st.append(adjacency_statistics)

    L_df = pd.DataFrame(laplacian_st)
    A_df = pd.DataFrame(adjacency_st)

    if inligprot:

        if expo:
            L_df.to_csv(name + 'L_stats_expo.csv', mode='a', header=False, index=False)
            A_df.to_csv(name + 'A_stats_expo.csv', mode='a', header=False, index=False)

        else:
            L_df.to_csv(name + 'L_stats_lorentz.csv', mode='a', header=False, index=False)
            A_df.to_csv(name + 'A_stats_lorentz.csv', mode='a', header=False, index=False)


    else:

        if expo:
            if os.path.isfile(name + "L_stats_expo.csv"):
                L_df.to_csv(name + 'L_stats_expo.csv', mode='a', header=False, index=False)
            else:
                with open(name + 'L_stats_expo.csv', "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(element_features)
                L_df.to_csv(name + 'L_stats_expo.csv', mode='a', header=False, index=False)


            if os.path.isfile(name + "A_stats_expo.csv"):
                A_df.to_csv(name + 'A_stats_expo.csv', mode='a', header=False, index=False)
            else:
                with open(name + 'A_stats_expo.csv', "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(element_features)
                A_df.to_csv(name + 'A_stats_expo.csv', mode='a', header=False, index=False)

        else:
            if os.path.isfile(name + "L_stats_lorentz.csv"):
                L_df.to_csv(name + 'L_stats_lorentz.csv', mode='a', header=False, index=False)
            else:
                with open(name + 'L_stats_lorentz.csv', "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(element_features)
                L_df.to_csv(name + 'L_stats_lorentz.csv', mode='a', header=False, index=False)

            if os.path.isfile(name + "A_stats_lorentz.csv"):
                A_df.to_csv(name + 'A_stats_lorentz.csv', mode='a', header=False, index=False)
            else:
                with open(name + 'A_stats_lorentz.csv', "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(element_features)
                A_df.to_csv(name + 'A_stats_lorentz.csv', mode='a', header=False, index=False)

    return 0


def main():

    start = time.time()

    path, ligand_name, at_types, selection_radius, geom_center, double_center, expo, lor, inligprot, name = parseArg()

    pose_store, prody_parsed = PDBParser(path=path)
    selected_protein, selected_ligand = binding_pocket_selection(pose_store=pose_store, p=prody_parsed, ligand_name=ligand_name, selection_radius=selection_radius, center=geom_center, double_center=double_center)

    if at_types:
        interactions, atom_types, ligand_atom_types, protein_atom_types = atomTypesDistanceCalc(binding_pocket=selected_protein, ligand=selected_ligand)
        final_weigths, atom_combinations = subgraphsWeights(atom_interactions=interactions, types=atom_types, expo=expo, lorentz=lor, inligprot=inligprot, ligand_atom_types=ligand_atom_types, protein_atom_types=protein_atom_types, at_types=at_types)
    else:
        interactions, elements, ligand_elements, protein_elements = elementsDistanceCalc(binding_pocket=selected_protein, ligand=selected_ligand, inligprot=inligprot)
        final_weigths, atom_combinations = subgraphsWeights(atom_interactions=interactions, types=elements, expo=expo, lorentz=lor, inligprot=inligprot, ligand_atom_types=ligand_elements, protein_atom_types=protein_elements, at_types=at_types)

    L_mat = laplacianMatrix(weights=final_weigths, atom_combinations=atom_combinations, at_types=at_types)
    A_mat = adjacencyMatrix(weights=final_weigths, atom_combinations=atom_combinations, at_types=at_types)
    LP = laplacianStats(matrices=L_mat, path=path)
    AD = adjacencyStats(matrices=A_mat, path=path)

    statsToCsv(laplacian_statistics=LP, adjacency_statistics=AD, expo=expo, lorentz=lor, inligprot=inligprot, at_types=at_types, name=name)

    end = time.time()

    print("\nRun time: {} seconds.".format(end - start))

    with open("run_times.txt", 'a') as rt:
        rt.write(str(end - start))
        rt.write('\n')

if __name__ == "__main__":
    main()
