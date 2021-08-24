import sys
import yaml
import time
import os

import argparse as ap

from keras.callbacks import ModelCheckpoint
import tensorflow as tf
from sklearn.metrics import r2_score

from parsing.parser import PDBParser, ligand_parse_write, read_graph_data

from graph.selection import binding_pocket_selection, ligand_atom_type_calc
from graph.distComp import elementsDistanceCalc, atomTypesDistanceCalc
from graph.graph_generator import atomSubgraphsWeights, elementSubgraphsWeights, laplacianMatrix, adjacencyMatrix
from graph.statistics import laplacianStats, adjacencyStats

from utils.utils import list_files, append_targets, statsToCsv, write_classification_report, write_regression_report
from utils.utils import write_regression_report_GBR, write_regression_report_XGBR, write_regression_report_LGBR, write_regression_report_MLPR, write_regression_report_ffnn

from preprocessing.preprocessing import dataset_preparation, data_splitting_regression, data_splitting_classification, Scaler, data_splitting_ffnn
from models.fit_graph import LearningModels


def parseyaml():

    with open(sys.argv[1]) as ctrl_file:
        params = yaml.load(ctrl_file, Loader=yaml.FullLoader)

    return params


def main():

    start = time.time()

    param_args = parseyaml()

    pipe = False
    try_all = True

    if 'cpus' not in param_args:
        cpus = 1

    if 'path_graph' not in param_args or pipe:

        pose_list = list_files(param_args['path'])

        for pose in pose_list:

            pose_store, prody_parsed = PDBParser(path=param_args['path'] + pose)
            selected_protein, selected_ligand = binding_pocket_selection(
                pose_store=pose_store, p=prody_parsed, ligand_name=param_args['ligand_name'], selection_radius=param_args['selection_radius'], center=param_args['center'])

            if param_args['nodes'] == 'atoms':
                ligand_path = ligand_parse_write(
                    path=param_args['path'] + pose, out=param_args['output'], lig_name=param_args['ligand_name'])
                selected_ligand_at = ligand_atom_type_calc(
                    ligand=selected_ligand, ligand_path=ligand_path)
                interactions, atom_types, ligand_atom_types, protein_atom_types = atomTypesDistanceCalc(
                    binding_pocket=selected_protein, ligand=selected_ligand_at)
                final_weigths, atom_combinations = atomSubgraphsWeights(atom_interactions=interactions, types=atom_types, decay_function=param_args['decay_function'],
                                                                        ligand_atom_types=ligand_atom_types, protein_atom_types=protein_atom_types)
            elif param_args['nodes'] == 'elements':
                interactions, elements, ligand_elements, protein_elements = elementsDistanceCalc(
                    binding_pocket=selected_protein, ligand=selected_ligand)
                final_weigths, atom_combinations = elementSubgraphsWeights(atom_interactions=interactions, types=elements, decay_function=param_args['decay_function'],
                                                                           ligand_atom_types=ligand_elements, protein_atom_types=protein_elements)

            L_mat = laplacianMatrix(weights=final_weigths,
                                    atom_combinations=atom_combinations)
            A_mat = adjacencyMatrix(weights=final_weigths,
                                    atom_combinations=atom_combinations)

            LP = laplacianStats(matrices=L_mat, pose=pose)
            AD = adjacencyStats(matrices=A_mat, pose=pose)

            statsToCsv(laplacian_statistics=LP, adjacency_statistics=AD,
                       decay_function=param_args['decay_function'], nodes=param_args['nodes'], name=param_args['run_name'], out=param_args['output'])

    if 'path_graph' in param_args or pipe:

        if not pipe:

            data = read_graph_data(path=param_args['path_graph'])

        elif pipe:

            data = append_targets(name=param_args['run_name'], out=param_args['output'],
                                  target=param_args['target'], decay_function=param_args['decay_function'])

        data = dataset_preparation(data=data)

        LM = LearningModels()

        if param_args['task'] == 'classification':
            train, activity_train, test, activity_test = data_splitting_classification(
                data=data, test_size=param_args['test_size'], seed=param_args['seed'])

            if 'algorithm' in param_args:

                if param_args['algorithm'] == 'GBC':

                    best_lr, best_train_acc, best_test_acc = LM.GBC_optimization(
                        train=train, activity_train=activity_train, test=test, activity_test=activity_test)
                    cf, report, MCC = LM.GBC(train=train, activity_train=activity_train,
                                             test=test, activity_test=activity_test, best_lr=best_lr)

                    classification_outputs = {'GBC_optimization': [best_lr, best_train_acc, best_test_acc],
                                              'GBC': [cf, report, MCC]}

                    write_classification_report_GBR(
                        classification_outputs=classification_outputs, out_file=param_args['output'] + 'classification_report.out')

                elif param_args['algorithm'] == 'XGBC':

                    XGB_accuracy, cf_XGB, report_XGB, MCC_XGB = LM.XGBoost(
                        train=train, activity_train=activity_train, test=test, activity_test=activity_test)

                    classification_outputs = {
                        'XGB': [XGB_accuracy, cf_XGB, report_XGB, MCC_XGB]}

                    write_classification_report_XGBR(
                        classification_outputs=classification_outputs, out_file=param_args['output'] + 'classification_report.out')

                elif param_args['algorithm'] == 'LGBC':

                    lgb_accuracy, cf_lgb, report_lgb, MCC_lgb = LM.LigthGB(
                        train=train, activity_train=activity_train, test=test, activity_test=activity_test)

                    classification_outputs = {
                        'LGBC': [lgb_accuracy, cf_lgb, report_lgb, MCC_lgb]}

                    write_classification_report_LGBR(
                        classification_outputs=classification_outputs, out_file=param_args['output'] + 'classification_report.out')

            else:

                best_lr, best_train_acc, best_test_acc = LM.GBC_optimization(
                    train=train, activity_train=activity_train, test=test, activity_test=activity_test)
                cf, report, MCC = LM.GBC(train=train, activity_train=activity_train,
                                         test=test, activity_test=activity_test, best_lr=best_lr)

                XGB_accuracy, cf_XGB, report_XGB, MCC_XGB = LM.XGBoost(
                    train=train, activity_train=activity_train, test=test, activity_test=activity_test)

                lgb_accuracy, cf_lgb, report_lgb, MCC_lgb = LM.LigthGB(
                    train=train, activity_train=activity_train, test=test, activity_test=activity_test)

                classification_outputs = {'GBC_optimization': [best_lr, best_train_acc, best_test_acc],
                                          'GBC': [cf, report, MCC], 'XGB': [XGB_accuracy, cf_XGB, report_XGB, MCC_XGB],
                                          'LGBC': [lgb_accuracy, cf_lgb, report_lgb, MCC_lgb]}

                write_classification_report(classification_outputs=classification_outputs,
                                            out_file=param_args['output'] + 'classification_report.out')

        elif param_args['task'] == 'regression':

            if param_args['algorithm'] == 'FFNN':

                S = Scaler()

                train, bindingEnergy_train, val, bindingEnergy_val, test, bindingEnergy_test,
                    ligands_test, ligandRMSD_test = data_splitting_ffnn(data=data,
                                                                        seed=param_args['seed'])

                train, val, test = S.min_max_scaler_ffnn(train=train, val=val, test=test)

                if param_args['pelePrep'] == 'profile':
                    NN_model = LM.FFNN_profile(param_args['learning_rate'])
                else:
                    NN_model = LM.FFNN_clustering(param_args['learning_rate'])

                checkpoint_name = 'weights.hdf5'
                checkpoint = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose = 1, save_best_only = True, mode ='auto')
                callbacks_list = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5), checkpoint]

                history = NN_model.fit(train, bindingEnergy_train, epochs=param_args['epochs'],
                                        seed=param_args['epochs'], batch_size=param_args['batch_size'],
                                        validation_data=(val, bindingEnergy_val), callbacks=callbacks_list)

                results_train = NN_model.evaluate(train, bindingEnergy_train)
                results_val = NN_model.evaluate(val, bindingEnergy_val)
                results_test = NN_model.evaluate(test, bindingEnergy_test)

                r2_train = r2_score(bindingEnergy_train, train_pred)
                r2_val = r2_score(bindingEnergy_val, val_pred)
                r2_test = r2_score(bindingEnergy_val, test_pred)

                results = [results_train.append(r2_train), results_val.append(r2_val),
                            results_test.append(r2_test)]

                write_regression_report_ffnn(results, param_args, param_args['output'])

                train_pred = NN_model.predict(train)
                val_pred = NN_model.predict(val)
                test_pred = NN_model.predict(test)

                train_plots(history, param_args['output'])

                if 'ligandRMSD' in list(data.columns):

                    target_plot(test_pred, ligand_test, bindingEnergy_test, ligandRMSD_test)

            else:

                if 'scaler' in param_args:

                    best_scaler = None
                    S = Scaler()

                    scaling_functions = [(S.standard_scaler, "standard"), (S.min_max_scaler, 'min_max'), (S.power_scaler, 'power'), (
                        S.quantile_scaler_uniform, 'uniform'), (S.quantile_scaler_gauss, 'quantile_Gauss'), (S.robust_scaler, 'robust'), (S.max_abs_scaler, 'max_abs')]

                    if param_args['scaler'] == 'search' and best_scaler is None:

                        scalers_performance = []

                        for scale in scaling_functions:

                            print("---------------------------------------")
                            print("\nApplying {} scaler\n".format(scale[1]))

                            try:
                                train, test = scale[0](train=train, test=test)

                            except ValueError:
                                print(
                                    "{} scaler cannot be applied to this data values.\n".format(scale[1]))
                                continue

                            pred, R2_test, MSE, MAE, epoch_pred = LM.GBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=0.25)

                            scalers_performance.append((R2_test, scale[0], scale[1]))

                        scalers_best_performance = sorted(scalers_performance, key=lambda tup: tup[0])

                        best_scaler = scalers_best_performance[0]

                    if (param_args['scaler'] == 'search' and best_scaler is not None) or (param_args['scaler'] != 'search'):

                        for scaler in scaling_functions:
                            if scaler[1] == param_args['scaler']:
                                scale = scaler
                                break

                        print("---------------------------------------")
                        print("\nApplying {} scaler\n".format(scale[1]))

                        try:
                            train, test = scale[0](train=train, test=test)

                        except ValueError:
                            print(
                                "{} scaler cannot be applied to this data values.\n".format(scale[1]))

                        if param_args['algorithm'] == 'GBR':

                            best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
                            pred, R2_test, MSE, MAE, epoch_pred = LM.GBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=best_lr)

                            regression_outputs = {'scaler': scale, 'GBR_optimization': [
                                best_lr, best_train_R2, best_test_R2], 'GBR': [R2_test, MSE, MAE, epoch_pred]}

                            write_regression_report_GBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'XGBR':

                            XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2 = LM.XGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            XGBR_pred, XGBR_R2_test, XGBR_MSE, XGBR_MAE, = LM.XGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=XGBR_best_lr, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'XGBR_optimization': [XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2],
                                                  'XGBR': [XGBR_R2_test, XGBR_MSE, XGBR_MAE]}

                            write_regression_report_XGBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'LGBR':

                            LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2 = LM.LGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE = LM.LGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=LGBR_best_lr, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'LGBR_optimization': [LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2],
                                                  'LGBR': [LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE]}

                            write_regression_report_LGBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'MLPR':

                            MLPR_pred, MLPR_R2_test, MLPR_MSE, MLRP_MAE, params = LM.MLPR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'MLPR': [
                                params, MLPR_pred, MLPR_R2_test, MLPR_MSE, MLPR_MAE]}

                            write_regression_report_MLPR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        else:

                            best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
                            pred, R2_test, MSE, MAE, epoch_pred = LM.GBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=best_lr)

                            XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2 = LM.XGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            XGBR_pred, XGBR_R2_test, XGBR_MSE, XGBR_MAE, = LM.XGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=XGBR_best_lr, cpus=param_args['cpus'])

                            LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2 = LM.LGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE = LM.LGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=LGBR_best_lr, cpus=param_args['cpus'])

                            MLPR_pred, MLPR_R2_test, MLPR_MSE, MLRP_MAE, params = LM.MLPR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'GBR_optimization': [best_lr, best_train_R2, best_test_R2], 'GBR': [R2_test, MSE, MAE, epoch_pred],
                                                  'XGBR_optimization': [XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2],
                                                  'XGBR': [XGBR_R2_test, XGBR_MSE, XGBR_MAE], 'LGBR_optimization': [LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2],
                                                  'LGBR': [LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE], 'MLPR': [params, MLPR_pred, MLPR_R2_test, MLPR_MSE, MLPR_MAE]}

                            write_regression_report(regression_outputs=regression_outputs,
                                                    out_file=param_args['output'] + 'regression_report.out')

                else:

                    print("\nFitting models without scaling...it may work but be careful!")

                    scale = None

                    if 'algorithm' in param_args:
                        if param_args['algorithm'] == 'GBR':

                            best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
                            pred, R2_test, MSE, MAE, epoch_pred = LM.GBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=best_lr)

                            regression_outputs = {'scaler': scale, 'GBR_optimization': [
                                best_lr, best_train_R2, best_test_R2], 'GBR': [R2_test, MSE, MAE, epoch_pred]}

                            write_regression_report_GBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'XGBR':

                            XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2 = LM.XGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            XGBR_pred, XGBR_R2_test, XGBR_MSE, XGBR_MAE, = LM.XGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=XGBR_best_lr, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'XGBR_optimization': [XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2],
                                                  'XGBR': [XGBR_R2_test, XGBR_MSE, XGBR_MAE]}

                            write_regression_report_XGBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'LGBR':

                            LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2 = LM.LGBR_optimization(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                            LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE = LM.LGBR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=LGBR_best_lr, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'LGBR_optimization': [LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2],
                                                  'LGBR': [LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE]}

                            write_regression_report_LGBR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                        elif param_args['algorithm'] == 'MLPR':

                            MLPR_pred, MLPR_R2_test, MLPR_MSE, MLRP_MAE, params = LM.MLPR(
                                train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])

                            regression_outputs = {'scaler': scale, 'MLPR': [
                                params, MLPR_pred, MLPR_R2_test, MLPR_MSE, MLPR_MAE]}

                            write_regression_report_MLPR(
                                regression_outputs=regression_outputs, out_file=param_args['output'] + 'regression_report.out')

                    else:

                        best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
                        pred, R2_test, MSE, MAE, epoch_pred = LM.GBR(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=best_lr)

                        XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2 = LM.XGBR_optimization(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                        XGBR_pred, XGBR_R2_test, XGBR_MSE, XGBR_MAE, = LM.XGBR(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=XGBR_best_lr, cpus=param_args['cpus'])

                        LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2 = LM.LGBR_optimization(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])
                        LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE = LM.LGBR(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, best_lr=LGBR_best_lr, cpus=param_args['cpus'])

                        MLPR_pred, MLPR_R2_test, MLPR_MSE, MLRP_MAE, params = LM.MLPR(
                            train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test, cpus=param_args['cpus'])

                        regression_outputs = {'scaler': scale, 'GBR_optimization': [best_lr, best_train_R2, best_test_R2], 'GBR': [R2_test, MSE, MAE, epoch_pred],
                                              'XGBR_optimization': [XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2],
                                              'XGBR': [XGBR_R2_test, XGBR_MSE, XGBR_MAE], 'LGBR_optimization': [LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2],
                                              'LGBR': [LGBR_pred, LGBR_R2_test, LGBR_MSE, LGBR_MAE], 'MLPR': [params, MLPR_pred, MLPR_R2_test, MLPR_MSE, MLPR_MAE]}

                        write_regression_report(regression_outputs=regression_outputs,
                                                out_file=param_args['output'] + 'regression_report.out')

    end = time.time()

    print("\nRun time: {} seconds.".format(end - start))

    with open(param_args['output'] + param_args['run_name'] + "_run_times.txt", 'a') as rt:
        rt.write(str(end - start))
        rt.write('\n')


if __name__ == "__main__":
    main()
