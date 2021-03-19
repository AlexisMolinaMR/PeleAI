import math
from math import sqrt


def elementsDistanceCalc(binding_pocket, ligand):
    '''
    Compute distances between selected elements and discard covalent interactions.
    '''

    elements = []
    at_count = 0

    atom_interactions = []

    used_elements = []

    atomic_radius_elements = {'C': 0.7, 'O': 0.6, 'N': 0.65, 'H': 0.25,
                              'Cl': 1, 'P': 1, 'S': 1, 'FE': 1.4, 'Br': 1.15, 'F': 0.5, 'I': 1.4}

    at_count = 0

    protein_elements = ['C', 'O', 'N', 'S']
    protein_elements.sort()
    ligand_elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
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
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j]
                                                [-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()
            except:
                if (x1 + y1 + z1 <= atomic_radius_elements[binding_pocket[i][-1]]**2) or (x1 + y1 + z1 <= atomic_radius_elements[ligand[i][-1]]**2):
                    print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                        protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                else:
                    distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j]
                                            [-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                    atom_interactions.append(distance_interaction)
                    distance_interaction = ()

    elements = None

    return atom_interactions, elements, ligand_elements, protein_elements


def atomTypesDistanceCalc(binding_pocket, ligand):

    atom_interactions = []

    atom_types = []

    used_types = []

    atomic_radius_types = {'NH1': 1.65, 'CP': 1.76, 'CH1E': 1.87, 'O': 0.6, 'OP': 1.4, 'CH0': 1.76, 'CH1S': 1.87, 'CH2E': 1.87, 'CH3E': 1.87, 'CR1E': 1.76, 'OH1': 1.4, 'OC': 1.4, 'OS': 1.4, 'CH2G': 1.87, 'CH2P': 1.87,
                           'NH1S': 1.65, 'NC2': 1.65, 'NH2': 1.65, 'CR1W': 1.76, 'CY2': 1.76, 'SC': 1.85, 'CF': 1.76, 'SM': 1.85, 'CY': 1.76, 'SM': 1.85, 'CY': 1.76, 'CW': 1.76, 'CRHH': 1.76, 'NH3': 1.5,
                           'CR1H': 1.76, 'C5': 1.76, 'C': 0.7, 'N': 1, 'NP': 1.65, 'C5W': 1.76, 'H': 0.25, 'Cl': 1, 'P': 1, 'S': 1, 'FE': 1.4, 'Br': 1.15, 'F': 0.5, 'I': 1.4}
    at_count = 0

    protein_atom_types = ['NH1', 'C', 'CP', 'CH1E', 'OP', 'O', 'CH0', 'CH1S', 'CH2E', 'CH3E', 'CR1E', 'OH1', 'OC', 'OS', 'CH2G', 'CH2P', 'NH1S', 'NC2', 'NH2', 'CR1W', 'CY2', 'SC', 'CF', 'SM', 'CY', 'CW', 'CRHH', 'NH3', 'CR1H', 'C5',
                          'N', 'NP', 'C5W', 'H', 'Cl', 'P', 'S', 'FE', 'Br', 'F', 'I']

    protein_atom_types.sort()
    ligand_atom_types = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
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
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (binding_pocket[i][-1], binding_pocket[i][1], ligand[j]
                                                [-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()
            except:
                if binding_pocket[i][1] == 'NH1':
                    if binding_pocket[i][0] == 'ARG':
                        if (x1 + y1 + z1 <= atomic_radius_types['NC2']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NC2', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()
                    else:
                        if (x1 + y1 + z1 <= atomic_radius_types['NH1']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NH1', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CA':
                    if binding_pocket[i][0] != 'GLY':
                        if (x1 + y1 + z1 <= atomic_radius_types['CH1E']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH1E', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    else:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH2G']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH2G', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1][0:2] == 'CZ':
                    if binding_pocket[i][0] == 'ARG':
                        if (x1 + y1 + z1 <= atomic_radius_types['CH0']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH0', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1][0:2] == 'CG':
                    if binding_pocket[i][0] in ['ASN', 'ASP']:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH0']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH0', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'LEU':
                        if (x1 + y1 + z1 <= atomic_radius_types['CH1S']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH1S', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'PHE':
                        if (x1 + y1 + z1 <= atomic_radius_types['CF']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CF', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'TYR':
                        if (x1 + y1 + z1 <= atomic_radius_types['CY']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CY', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'HIS':
                        if (x1 + y1 + z1 <= atomic_radius_types['C5']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'C5', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'TRP':
                        if (x1 + y1 + z1 <= atomic_radius_types['C5W']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'C5W', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    else:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH2P']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH2P', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CD':
                    if binding_pocket[i][0] in ['GLN', 'GLU']:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH0']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH0', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    else:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH2P']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH2P', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CB':
                    if binding_pocket[i][0] in ['ILE', 'THR', 'VAL']:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH1S']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH1S', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'PRO':
                        if (x1 + y1 + z1 <= atomic_radius_types['CH2P']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH2P', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CH2':
                    if binding_pocket[i][0] not in ['GLY', 'PRO']:
                        if (x1 + y1 + z1 <= atomic_radius_types['CH2E']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CH2E', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()
                    else:
                        if (x1 + y1 + z1 <= atomic_radius_types['CR1W']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CR1W', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CH3':
                    if (x1 + y1 + z1 <= atomic_radius_types['CH3E']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'CH3E', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] == 'CH':
                    if (x1 + y1 + z1 <= atomic_radius_types['CR1E']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'CR1E', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] in ['OG', 'OG1', 'OH']:
                    if binding_pocket[i][0] not in ['SER', 'THR', 'TYR']:
                        if (x1 + y1 + z1 <= atomic_radius_types['OH1']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'OH1', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] in ['OD2', 'OE2']:
                    if (x1 + y1 + z1 <= atomic_radius_types['OC']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'OC', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] in ['OD1', 'OE1']:
                    if binding_pocket[0] in ['ASP', 'GLU']:
                        if (x1 + y1 + z1 <= atomic_radius_types['OC']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'OC', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    if binding_pocket[i][0] in ['ASN', 'GLN']:
                        if (x1 + y1 + z1 <= atomic_radius_types['OS']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'OS', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'NE1':
                    if (x1 + y1 + z1 <= atomic_radius_types['NH1S']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'NH1S', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] in ['NE', 'ND1']:
                    if binding_pocket[i][0] in ['ARG', 'HIS']:
                        if (x1 + y1 + z1 <= atomic_radius_types['NH1S']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NH1S', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'NH2':
                    if (x1 + y1 + z1 <= atomic_radius_types['NC2']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'NC2', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] == 'NH1':
                    if binding_pocket[i][0] == 'ARG':
                        if (x1 + y1 + z1 <= atomic_radius_types['NC2']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NC2', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] in ['ND2', 'NE2']:
                    if binding_pocket[i][0] in ['ASN', 'GLN']:
                        if (x1 + y1 + z1 <= atomic_radius_types['NH2']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NH2', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CZ2':
                    if binding_pocket[i][0] == 'TRP':
                        if (x1 + y1 + z1 <= atomic_radius_types['CR1W']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CR1W', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CZ':
                    if binding_pocket[i][0] == 'TYR':
                        if (x1 + y1 + z1 <= atomic_radius_types['CY2']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CY2', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'S':
                    if binding_pocket[i][0] == 'CYS':
                        if (x1 + y1 + z1 <= atomic_radius_types['SC']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'SC', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'MET':
                        if (x1 + y1 + z1 <= atomic_radius_types['SM']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'SM', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CE2':
                    if (x1 + y1 + z1 <= atomic_radius_types['CW']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'CW', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] == 'CD2':
                    if binding_pocket[i][0] == 'TRP':
                        if (x1 + y1 + z1 <= atomic_radius_types['CW']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CW', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                    elif binding_pocket[i][0] == 'HIS':
                        if (x1 + y1 + z1 <= atomic_radius_types['CR1H']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CR1H', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'CE1':
                    if binding_pocket[i][0] == 'HIS':
                        if (x1 + y1 + z1 <= atomic_radius_types['CRHH']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'CRHH', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'NZ':
                    if binding_pocket[i][0] == 'LYS':
                        if (x1 + y1 + z1 <= atomic_radius_types['NH3']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NH3', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'N':
                    if binding_pocket[i][0] == 'PRO':
                        if (x1 + y1 + z1 <= atomic_radius_types['NP']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                            print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                                protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                        else:
                            distance_interaction = (
                                binding_pocket[i][-1], 'NP', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                            atom_interactions.append(distance_interaction)
                            distance_interaction = ()

                elif binding_pocket[i][1] == 'C':
                    if (x1 + y1 + z1 <= atomic_radius_types['CP']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'CP', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1] == 'O':
                    if (x1 + y1 + z1 <= atomic_radius_types['OP']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'OP', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

                elif binding_pocket[i][1][0] == 'H':
                    if (x1 + y1 + z1 <= atomic_radius_types['H']**2) or (x1 + y1 + z1 <= atomic_radius_types[ligand[i][-1]]**2):
                        print("Distance of between atoms {} ({}) and {} ({}) from residues {} and {} -> Covalent interaction or in same residue".format(
                            protein_ligand[i][-1], protein_ligand[i][1], protein_ligand[j][-1], protein_ligand[j][1], protein_ligand[i][0], protein_ligand[j][0], math.sqrt(x1 + y1 + z1)))

                    else:
                        distance_interaction = (
                            binding_pocket[i][-1], 'OP', ligand[j][-1], ligand[j][1], binding_pocket[i][0], ligand[j][0], math.sqrt(x1 + y1 + z1))
                        atom_interactions.append(distance_interaction)
                        distance_interaction = ()

    for i in sorted(protein_atom_types):
        for j in sorted(ligand_atom_types):
            atom_types.append(j + '-' + i)

    atom_types.sort()

    return atom_interactions, atom_types, ligand_atom_types, protein_atom_types
