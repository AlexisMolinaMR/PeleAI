import csv
import os

import pandas as pd

from os import listdir


def list_files(path):
    pose_list = [f for f in listdir(path) if f.endswith('.' + 'pdb')]
    return pose_list


def append_targets(name, out, target, decay_function):

    features = pd.read_csv(out + name + '_' + 'L_stats_' + decay_function + '.csv')
    targets = pd.read_csv(target)

    features_targets = pd.merge(features, targets, on='ligand')

    return features_targets


def statsToCsv(laplacian_statistics, adjacency_statistics, decay_function, nodes, name, out):
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
        if os.path.isfile(out + name + '_' + "L_stats_expo.csv"):
            L_df.to_csv(out + name + '_' + 'L_stats_expo.csv', mode='a', header=False, index=False)
        else:
            with open(out + name + '_' + 'L_stats_expo.csv', "w") as f:
                writer = csv.writer(f)
                writer.writerow(element_features)
            L_df.to_csv(out + name + '_' + 'L_stats_expo.csv', mode='a', header=False, index=False)

        if os.path.isfile(out + name + '_' + "A_stats_expo.csv"):
            A_df.to_csv(out + name + '_' + 'A_stats_expo.csv', mode='a', header=False, index=False)
        else:
            with open(out + name + '_' + 'A_stats_expo.csv', "w") as f:
                writer = csv.writer(f)
                writer.writerow(element_features)
            A_df.to_csv(out + name + '_' + 'A_stats_expo.csv', mode='a', header=False, index=False)

    elif decay_function == 'lorentz':
        if os.path.isfile(out + name + '_' + "L_stats_lorentz.csv"):
            L_df.to_csv(out + name + '_' + 'L_stats_lorentz.csv',
                        mode='a', header=False, index=False)
        else:
            with open(out + name + '_' + 'L_stats_lorentz.csv', "w") as f:
                writer = csv.writer(f)
                writer.writerow(element_features)
            L_df.to_csv(out + name + '_' + 'L_stats_lorentz.csv',
                        mode='a', header=False, index=False)

        if os.path.isfile(out + name + '_' + "A_stats_lorentz.csv"):
            A_df.to_csv(out + name + '_' + 'A_stats_lorentz.csv',
                        mode='a', header=False, index=False)
        else:
            with open(out + name + '_' + 'A_stats_lorentz.csv', "w") as f:
                writer = csv.writer(f)
                writer.writerow(element_features)
            A_df.to_csv(out + name + '_' + 'A_stats_lorentz.csv',
                        mode='a', header=False, index=False)

    return 0


def write_classification_report_GBR(classification_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print("Gradient Boosting Classifier\n", file=out)
        print("Best learning rate: {}".format(
            classification_outputs['GBC_optimization'][0]), file=out)
        print("Best train accuracy: {}".format(
            classification_outputs['GBC_optimization'][1]), file=out)
        print("Best test accuracy: {}\n".format(
            classification_outputs['GBC_optimization'][2]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['GBC'][0], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['GBC'][1], file=out)
        print("\nMCC: {}".format(classification_outputs['GBC'][2]), file=out)

    return 0


def write_classification_report_XGBR(classification_outputs, out_file):

    with open(out_file, 'a') as out:

        print("\n#####################################", file=out)
        print("XGBoost Classifier\n", file=out)
        print("XGB accuracy: {}\n".format(classification_outputs['XGB'][0]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['XGB'][1], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['XGB'][2], file=out)
        print("\nMCC: {}".format(classification_outputs['XGB'][3]), file=out)

    return 0


def write_classification_report_LGBR(classification_outputs, out_file):

    with open(out_file, 'a') as out:

        print("\n#####################################", file=out)
        print("LightGBoost Classifier\n", file=out)
        print("LGBMC accuracy: {}\n".format(classification_outputs['LGBC'][0]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['LGBC'][1], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['LGBC'][2], file=out)
        print("\nMCC: {}".format(classification_outputs['LGBC'][2]), file=out)

    return 0


def write_classification_report(classification_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print("Gradient Boosting Classifier\n", file=out)
        print("Best learning rate: {}".format(
            classification_outputs['GBC_optimization'][0]), file=out)
        print("Best train accuracy: {}".format(
            classification_outputs['GBC_optimization'][1]), file=out)
        print("Best test accuracy: {}\n".format(
            classification_outputs['GBC_optimization'][2]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['GBC'][0], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['GBC'][1], file=out)
        print("\nMCC: {}".format(classification_outputs['GBC'][2]), file=out)
        print("\n#####################################", file=out)
        print("XGBoost Classifier\n", file=out)
        print("XGB accuracy: {}\n".format(classification_outputs['XGB'][0]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['XGB'][1], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['XGB'][2], file=out)
        print("\nMCC: {}".format(classification_outputs['XGB'][3]), file=out)
        print("\n#####################################", file=out)
        print("LightGBoost Classifier\n", file=out)
        print("LGBMC accuracy: {}\n".format(classification_outputs['LGBC'][0]), file=out)
        print("Confusion matrix", file=out)
        print(classification_outputs['LGBC'][1], file=out)
        print("\nClassification report:", file=out)
        print(classification_outputs['LGBC'][2], file=out)
        print("\nMCC: {}".format(classification_outputs['LGBC'][2]), file=out)

    return 0


def write_regression_report_GBR(regression_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print('Scaler: {}'.format(regression_outputs['scaler']), file=out)
        print("Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(regression_outputs['GBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['GBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['GBR'][1]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['GBR'][2]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['GBR'][3]), file=out)

    return 0


def write_regression_report_XGBR(regression_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print("Xtreme Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(
            regression_outputs['XGBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['XGBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['XGBR'][0]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['XGBR'][1]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['XGBR'][2]), file=out)

    return 0


def write_regression_report_LGBR(regression_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print("Light Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(
            regression_outputs['LGBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['LGBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['LGBR'][1]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['LGBR'][2]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['XGBR'][3]), file=out)

    return 0


def write_regression_report_MLPR(regression_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print("Multi-layer Perceptron Regressor\n", file=out)
        print("Best parameters: {}".format(regression_outputs['MLPR'][0]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['MLPR'][2]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['MLPR'][3]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['MLPR'][4]), file=out)

    return 0


def write_regression_report(regression_outputs, out_file):

    with open(out_file, 'a') as out:

        print("#####################################", file=out)
        print('Scaler: {}'.format(regression_outputs['scaler']), file=out)
        print("Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(regression_outputs['GBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['GBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['GBR'][1]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['GBR'][2]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['GBR'][3]), file=out)
        print("#####################################", file=out)
        print("Xtreme Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(
            regression_outputs['XGBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['XGBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['XGBR'][0]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['XGBR'][1]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['XGBR'][2]), file=out)
        print("#####################################", file=out)
        print("Light Gradient Boosting Regressor\n", file=out)
        print("Best learning rate: {}".format(
            regression_outputs['LGBR_optimization'][0]), file=out)
        print("Best train R2: {}".format(regression_outputs['LGBR_optimization'][1]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['LGBR'][1]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['LGBR'][2]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['XGBR'][3]), file=out)
        print("#####################################", file=out)
        print("Multi-layer Perceptron Regressor\n", file=out)
        print("Best parameters: {}".format(regression_outputs['MLPR'][0]), file=out)
        print("Best test R2: {}\n".format(regression_outputs['MLPR'][2]), file=out)
        print("Mean squared error: {}\n".format(regression_outputs['MLPR'][3]), file=out)
        print("Mean absolute error: {}\n".format(regression_outputs['MLPR'][4]), file=out)

    return 0
