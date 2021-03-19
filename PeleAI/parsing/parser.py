import pandas as pd

from prody import parsePDB


def PDBParser(path):
    '''
    This function will parse the different docking poses
    in PDB format.
    '''

    pose_store = {}

    with open(path, 'r') as pdb:
        for lines in pdb:
            line = lines.strip()
            if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20].strip() in ["MSE", "KCX"]):
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


def read_graph_data(path):
    '''
    Read csv file with features extracted for each pose.
    '''

    d = pd.read_csv(path)
    data = d.copy()

    return data
