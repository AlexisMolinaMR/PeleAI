import os
import math

import numpy as np
import pandas as pd


def elementSubgraphsWeights(atom_interactions, types, decay_function, ligand_atom_types, protein_atom_types):
    '''
    Compute subgraph weights using either exponential or Lorentz decay functions.
    '''

    weights = []
    count = 0
    atom_atom = []
    atom_rev = []
    uni_inter = []

    final_weigths = []

    if decay_function == 'expo':
        for inter in atom_interactions:
            type_weight = (np.exp(-abs(inter[-1])/3), "{}-{}".format(inter[0], inter[2]))
            weights.append(type_weight)
            type_weight = ()

    elif decay_function == 'lorentz':
        for inter in atom_interactions:
            type_weight = ((1/1 + (abs(inter[-1])/3)**5), "{}-{}".format(inter[0], inter[2]))
            weights.append(type_weight)
            type_weight = ()

    atom_atom.sort()

    return final_weigths, atom_atom


def atomSubgraphsWeights(atom_interactions, types, decay_function, ligand_atom_types, protein_atom_types):
    '''
    Compute subgraph weights using either exponential or Lorentz decay functions.
    '''

    weights = []
    count = 0
    atom_atom = []
    atom_rev = []
    uni_inter = []

    final_weigths = []

    if decay_function == 'expo':
        for inter in atom_interactions:
            type_weight = (np.exp(-abs(inter[-1])/3), "{}-{}".format(inter[2], inter[1]))
            weights.append(type_weight)
            type_weight = ()

    elif decay_function == 'lorentz':
        for inter in atom_interactions:
            type_weight = ((1/1 + (abs(inter[-1])/3)**5), "{}-{}".format(inter[2], inter[1]))
            weights.append(type_weight)
            type_weight = ()

    for i in range(len(ligand_atom_types)):
        for j in range(len(protein_atom_types)):
            atom_atom.append(("{}-{}".format(ligand_atom_types[i], protein_atom_types[j])))
            if ligand_atom_types[i] != protein_atom_types[j] and "{}-{}".format(ligand_atom_types[i], protein_atom_types[j]) not in atom_rev:
                atom_rev.append(("{}-{}".format(protein_atom_types[j], ligand_atom_types[i])))

    weights = list(set(weights))

    final_weigths = weights

    atom_atom.sort()

    return final_weigths, atom_atom


def laplacianMatrix(weights, atom_combinations):
    '''
    Compute Laplacian matrices.
    '''

    laplacian_matrices = []

    for element in atom_combinations:
        element_count = 0
        for weight in weights:
            if weight[-1] == element:
                element_count += 1
        print("\nLaplacian matrix for subgraph G({}). Matrix size {}x{}:\n".format(
            element, element_count, element_count))
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


def adjacencyMatrix(weights, atom_combinations):
    '''
    Compute adjacency matrices.
    '''

    adjacency_matrices = []

    for element in atom_combinations:
        element_count = 0
        for weight in weights:
            if weight[-1] == element:
                element_count += 1
        print("\nAdjacency matrix for subgraph G({}). Matrix size {}x{}:\n".format(
            element, element_count, element_count))
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


def laplacianStats(matrices, path):
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
            total_sec_pow = 0
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


def statsToCsv(laplacian_statistics, adjacency_statistics, decay_function, nodes, name):
    '''
    Write previously computed statistics to a csv file.
    '''

    if nodes == 'atoms':

        element_element = ['Br-Br', 'Br-C', 'Br-C5', 'Br-C5W', 'Br-CF', 'Br-CH0', 'Br-CH1E', 'Br-CH1S', 'Br-CH2E', 'Br-CH2G', 'Br-CH2P', 'Br-CH3E', 'Br-CP', 'Br-CR1E', 'Br-CR1H', 'Br-CR1W', 'Br-CRHH', 'Br-CW', 'Br-CY', 'Br-CY2', 'Br-Cl', 'Br-F', 'Br-FE', 'Br-H', 'Br-I', 'Br-N', 'Br-NC2', 'Br-NH1', 'Br-NH1S', 'Br-NH2', 'Br-NH3', 'Br-NP', 'Br-O', 'Br-OC', 'Br-OH1', 'Br-OP', 'Br-OS', 'Br-P', 'Br-S', 'Br-SC', 'Br-SM', 'C-Br', 'C-C', 'C-C5', 'C-C5W', 'C-CF', 'C-CH0', 'C-CH1E', 'C-CH1S', 'C-CH2E', 'C-CH2G', 'C-CH2P', 'C-CH3E', 'C-CP', 'C-CR1E', 'C-CR1H', 'C-CR1W', 'C-CRHH', 'C-CW', 'C-CY', 'C-CY2', 'C-Cl', 'C-F', 'C-FE', 'C-H', 'C-I', 'C-N', 'C-NC2', 'C-NH1', 'C-NH1S', 'C-NH2', 'C-NH3', 'C-NP', 'C-O', 'C-OC', 'C-OH1', 'C-OP', 'C-OS', 'C-P', 'C-S', 'C-SC', 'C-SM', 'Cl-Br', 'Cl-C', 'Cl-C5', 'Cl-C5W', 'Cl-CF', 'Cl-CH0', 'Cl-CH1E', 'Cl-CH1S', 'Cl-CH2E', 'Cl-CH2G', 'Cl-CH2P', 'Cl-CH3E', 'Cl-CP', 'Cl-CR1E', 'Cl-CR1H', 'Cl-CR1W', 'Cl-CRHH', 'Cl-CW', 'Cl-CY', 'Cl-CY2', 'Cl-Cl', 'Cl-F', 'Cl-FE', 'Cl-H', 'Cl-I', 'Cl-N', 'Cl-NC2', 'Cl-NH1', 'Cl-NH1S', 'Cl-NH2', 'Cl-NH3', 'Cl-NP', 'Cl-O', 'Cl-OC', 'Cl-OH1', 'Cl-OP', 'Cl-OS', 'Cl-P', 'Cl-S', 'Cl-SC', 'Cl-SM', 'F-Br', 'F-C', 'F-C5', 'F-C5W', 'F-CF', 'F-CH0', 'F-CH1E', 'F-CH1S', 'F-CH2E', 'F-CH2G', 'F-CH2P', 'F-CH3E', 'F-CP', 'F-CR1E', 'F-CR1H', 'F-CR1W', 'F-CRHH', 'F-CW', 'F-CY', 'F-CY2', 'F-Cl', 'F-F', 'F-FE', 'F-H', 'F-I', 'F-N', 'F-NC2', 'F-NH1', 'F-NH1S', 'F-NH2', 'F-NH3', 'F-NP', 'F-O', 'F-OC', 'F-OH1', 'F-OP', 'F-OS', 'F-P', 'F-S', 'F-SC', 'F-SM', 'H-Br', 'H-C', 'H-C5', 'H-C5W', 'H-CF', 'H-CH0', 'H-CH1E', 'H-CH1S', 'H-CH2E', 'H-CH2G', 'H-CH2P', 'H-CH3E', 'H-CP', 'H-CR1E', 'H-CR1H', 'H-CR1W', 'H-CRHH', 'H-CW', 'H-CY', 'H-CY2', 'H-Cl', 'H-F', 'H-FE', 'H-H', 'H-I', 'H-N', 'H-NC2', 'H-NH1', 'H-NH1S', 'H-NH2', 'H-NH3', 'H-NP', 'H-O', 'H-OC', 'H-OH1',
                           'H-OP', 'H-OS', 'H-P', 'H-S', 'H-SC', 'H-SM', 'I-Br', 'I-C', 'I-C5', 'I-C5W', 'I-CF', 'I-CH0', 'I-CH1E', 'I-CH1S', 'I-CH2E', 'I-CH2G', 'I-CH2P', 'I-CH3E', 'I-CP', 'I-CR1E', 'I-CR1H', 'I-CR1W', 'I-CRHH', 'I-CW', 'I-CY', 'I-CY2', 'I-Cl', 'I-F', 'I-FE', 'I-H', 'I-I', 'I-N', 'I-NC2', 'I-NH1', 'I-NH1S', 'I-NH2', 'I-NH3', 'I-NP', 'I-O', 'I-OC', 'I-OH1', 'I-OP', 'I-OS', 'I-P', 'I-S', 'I-SC', 'I-SM', 'N-Br', 'N-C', 'N-C5', 'N-C5W', 'N-CF', 'N-CH0', 'N-CH1E', 'N-CH1S', 'N-CH2E', 'N-CH2G', 'N-CH2P', 'N-CH3E', 'N-CP', 'N-CR1E', 'N-CR1H', 'N-CR1W', 'N-CRHH', 'N-CW', 'N-CY', 'N-CY2', 'N-Cl', 'N-F', 'N-FE', 'N-H', 'N-I', 'N-N', 'N-NC2', 'N-NH1', 'N-NH1S', 'N-NH2', 'N-NH3', 'N-NP', 'N-O', 'N-OC', 'N-OH1', 'N-OP', 'N-OS', 'N-P', 'N-S', 'N-SC', 'N-SM', 'O-Br', 'O-C', 'O-C5', 'O-C5W', 'O-CF', 'O-CH0', 'O-CH1E', 'O-CH1S', 'O-CH2E', 'O-CH2G', 'O-CH2P', 'O-CH3E', 'O-CP', 'O-CR1E', 'O-CR1H', 'O-CR1W', 'O-CRHH', 'O-CW', 'O-CY', 'O-CY2', 'O-Cl', 'O-F', 'O-FE', 'O-H', 'O-I', 'O-N', 'O-NC2', 'O-NH1', 'O-NH1S', 'O-NH2', 'O-NH3', 'O-NP', 'O-O', 'O-OC', 'O-OH1', 'O-OP', 'O-OS', 'O-P', 'O-S', 'O-SC', 'O-SM', 'P-Br', 'P-C', 'P-C5', 'P-C5W', 'P-CF', 'P-CH0', 'P-CH1E', 'P-CH1S', 'P-CH2E', 'P-CH2G', 'P-CH2P', 'P-CH3E', 'P-CP', 'P-CR1E', 'P-CR1H', 'P-CR1W', 'P-CRHH', 'P-CW', 'P-CY', 'P-CY2', 'P-Cl', 'P-F', 'P-FE', 'P-H', 'P-I', 'P-N', 'P-NC2', 'P-NH1', 'P-NH1S', 'P-NH2', 'P-NH3', 'P-NP', 'P-O', 'P-OC', 'P-OH1', 'P-OP', 'P-OS', 'P-P', 'P-S', 'P-SC', 'P-SM', 'S-Br', 'S-C', 'S-C5', 'S-C5W', 'S-CF', 'S-CH0', 'S-CH1E', 'S-CH1S', 'S-CH2E', 'S-CH2G', 'S-CH2P', 'S-CH3E', 'S-CP', 'S-CR1E', 'S-CR1H', 'S-CR1W', 'S-CRHH', 'S-CW', 'S-CY', 'S-CY2', 'S-Cl', 'S-F', 'S-FE', 'S-H', 'S-I', 'S-N', 'S-NC2', 'S-NH1', 'S-NH1S', 'S-NH2', 'S-NH3', 'S-NP', 'S-O', 'S-OC', 'S-OH1', 'S-OP', 'S-OS', 'S-P', 'S-S', 'S-SC', 'S-SM']

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

    elif nodes == 'elements':

        element_element = ['Br-C', 'Br-N', 'Br-O', 'Br-S', 'C-C', 'C-N', 'C-O', 'C-S', 'Cl-C', 'Cl-N', 'Cl-O', 'Cl-S', 'F-C', 'F-N', 'F-O', 'F-S', 'H-C', 'H-N',
                           'H-O', 'H-S', 'I-C', 'I-N', 'I-O', 'I-S', 'N-C', 'N-N', 'N-O', 'N-S', 'O-C', 'O-N', 'O-O', 'O-S', 'P-C', 'P-N', 'P-O', 'P-S', 'S-C', 'S-N', 'S-O', 'S-S']
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

    if decay_function == 'expo':
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

    elif decay_function == 'lorentz':
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
    selected_protein, selected_ligand = binding_pocket_selection(
        pose_store=pose_store, p=prody_parsed, ligand_name=ligand_name, selection_radius=selection_radius, center=geom_center, double_center=double_center)

    if at_types:
        interactions, atom_types, ligand_atom_types, protein_atom_types = atomTypesDistanceCalc(
            binding_pocket=selected_protein, ligand=selected_ligand)
        final_weigths, atom_combinations = subgraphsWeights(atom_interactions=interactions, types=atom_types, expo=expo, lorentz=lor,
                                                            inligprot=inligprot, ligand_atom_types=ligand_atom_types, protein_atom_types=protein_atom_types, at_types=at_types)
    else:
        interactions, elements, ligand_elements, protein_elements = elementsDistanceCalc(
            binding_pocket=selected_protein, ligand=selected_ligand, inligprot=inligprot)
        final_weigths, atom_combinations = subgraphsWeights(atom_interactions=interactions, types=elements, expo=expo, lorentz=lor,
                                                            inligprot=inligprot, ligand_atom_types=ligand_elements, protein_atom_types=protein_elements, at_types=at_types)

    L_mat = laplacianMatrix(weights=final_weigths,
                            atom_combinations=atom_combinations, at_types=at_types)
    A_mat = adjacencyMatrix(weights=final_weigths,
                            atom_combinations=atom_combinations, at_types=at_types)
    LP = laplacianStats(matrices=L_mat, path=path)
    AD = adjacencyStats(matrices=A_mat, path=path)

    statsToCsv(laplacian_statistics=LP, adjacency_statistics=AD, expo=expo,
               lorentz=lor, inligprot=inligprot, at_types=at_types, name=name)

    end = time.time()

    print("\nRun time: {} seconds.".format(end - start))

    with open("run_times.txt", 'a') as rt:
        rt.write(str(end - start))
        rt.write('\n')


if __name__ == "__main__":
    main()
