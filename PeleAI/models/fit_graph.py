from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, matthews_corrcoef, mean_squared_error
from sklearn.model_selection import StratifiedKFold, train_test_split, GridSearchCV, RandomizedSearchCV
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler, RobustScaler, Normalizer, QuantileTransformer, PowerTransformer
from sklearn import feature_selection
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier, XGBRegressor
from lightgbm import LGBMClassifier, LGBMRegressor

from math import sqrt

import pandas as pd
import numpy as np
import re

import time
import argparse as ap
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import warnings
warnings.filterwarnings("ignore")

global scaling_functions
scaling_functions = []

def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--path", required=True, type=str, help="path to csv file containing PeleAI data")
    parser.add_argument("-t", "--test_size", required=True, type=float, help="size of test set")
    parser.add_argument("-s", "--seed", required=True, type=int, help="random seed for train/test split")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', "--classification", action='store_true', help="fit into a classification algorithm")
    group.add_argument('-r', "--regression", action='store_true', help="fit into a regression algorithm")
    parser.add_argument("-d", "--data_scaling", action='store_true', help="try data scaling")

    args = parser.parse_args()

    path = args.path
    test_size = args.test_size
    seed = args.seed
    classification = args.classification
    regression = args.regression
    data_scaling = args.data_scaling

    return path, test_size, seed, classification, regression, data_scaling



def read_graph_data(path):
    '''
    Read csv file with features extracted for each pose.
    '''

    d = pd.read_csv(path)
    data = d.copy()

    return data


def dataset_preparation(data):
    '''
    Prepare dataframe for further data preparation.
    '''

    data.reset_index(drop=True, inplace= True)
    data.fillna(0, inplace=True)

    return data



def data_splitting_classification(data, test_size, seed):
    '''
    Split the date in train and test sets if classification problem.
    '''

    train, test = train_test_split(data, test_size=test_size, random_state=seed)

    activity_train = train['activity']
    activity_test = test['activity']

    train.drop(['ligand','activity'],axis=1, inplace=True)
    test.drop(['ligand','activity'], axis=1, inplace=True)

    return train, activity_train, test, activity_test

def data_splitting_regression(data, test_size, seed):
    '''
    Split the date in train and test sets if regression problem.
    '''

    train, test = train_test_split(data, test_size=test_size, random_state=seed)

    bindingEnergy_train = train['bindingEnergy']
    bindingEnergy_test = test['bindingEnergy']

    ligands_test = test['ligand']

    train.drop(['ligand','bindingEnergy'],axis=1, inplace=True)
    test.drop(['ligand','bindingEnergy'], axis=1, inplace=True)

    return train, bindingEnergy_train, test, bindingEnergy_test, ligands_test

class Scaler():
    '''
    Class containing all methods for scaling train and test data.
    '''

    @staticmethod
    def standard_scaler(train, test):
        '''
        Standardize features by removing the mean and scaling to unit variance.
        '''
        scaler = StandardScaler()

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def min_max_scaler(train, test):
        '''
        Transform features by scaling each feature to a given range. Set on default (0,1).
        '''

        scaler = MinMaxScaler()

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def max_abs_scaler(train, test):
        '''
        Scale each feature by its maximum absolute value.
        '''

        scaler = MaxAbsScaler()

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def robust_scaler(train, test):
        '''
        Scaler removes the median and scales the data according
        to the quantile range (defaults to IQR: Interquartile Range).
        '''

        scaler = RobustScaler()

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def normal_scaler(train, test):
        '''
        Normalize samples individually to unit norm.
        '''

        scaler = Normalizer()

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def quantile_scaler_gauss(train, test):
        '''
        Transforms the features to follow a normal distribution.
        '''

        scaler = QuantileTransformer(output_distribution='normal')

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def quantile_scaler_uniform(train, test):
        '''
        Transforms the features to follow a uniform distribution.
        '''

        scaler = QuantileTransformer(output_distribution='uniform')

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test

    @staticmethod
    def power_scaler(train, test):
        '''
        Apply a power transform featurewise to make data more Gaussian-like.
        '''

        scaler = PowerTransformer(method='yeo-johnson')

        train = scaler.fit_transform(train)
        test = scaler.fit_transform(test)

        return train, test



class LearningModels():
    '''
    This class contains all learning algorithms used for fitting the data.
    '''
    @staticmethod
    def GBC_optimization(train, activity_train, test, activity_test):
        '''
        Parameter optimization for the Gradient Boosting Classifier.
        '''

        lr_list = [0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
        accuracies_train = []
        accuracies_test = []
        best_test_acc = None
        best_train_acc = None
        best_lr = None
        for learning_rate in lr_list:
            gb_clf = GradientBoostingClassifier(n_estimators=20, learning_rate=learning_rate, max_features=2, max_depth=2, random_state=0)
            gb_clf.fit(train, activity_train)

            accuracies_train.append(gb_clf.score(train, activity_train))
            accuracies_test.append(gb_clf.score(test, activity_test))

        best_test_acc = max(accuracies_test)
        best_train_acc = accuracies_train[accuracies_test.index(best_test_acc)]
        best_lr = lr_list[accuracies_test.index(best_test_acc)]

        return best_lr, best_train_acc, best_test_acc

    @staticmethod
    def GBC(best_lr):
        '''
        GradientBoostingClassifier algorithm.
        '''

        gb_clf2 = GradientBoostingClassifier(n_estimators=20, learning_rate=best_lr, max_features=2, max_depth=2, random_state=0)
        gb_clf2.fit(train, activity_train)
        predictions = gb_clf2.predict(test)

        cf = confusion_matrix(activity_test, predictions)
        report = classification_report(activity_test, predictions)
        MCC = matthews_corrcoef(activity_test, predictions)

        return cf, report, MCC

    @staticmethod
    def XGBoost(train, activity_train, test, activity_test):
        '''
        XGBClassifier algorithm.
        '''

        xgb_clf = XGBClassifier()

        xgb_clf.fit(train, activity_train)

        y_pred = xgb_clf.predict(test)
        predictions = [round(value) for value in y_pred]

        XGB_accuracy = accuracy_score(activity_test, predictions)

        cf_XGB = confusion_matrix(activity_test, predictions)
        report_XGB = classification_report(activity_test, predictions)
        MCC_XGB = matthews_corrcoef(activity_test, predictions)

        return XGB_accuracy, cf_XGB, report_XGB, MCC_XGB

    @staticmethod
    def LigthGB(train, activity_train, test, activity_test):
        '''
        LigthGB algorithm.
        '''

        train_light = train.rename(columns = lambda x:re.sub('[^A-Za-z0-9_]+', '', x))
        test_light = test.rename(columns = lambda x:re.sub('[^A-Za-z0-9_]+', '', x))

        lgb_clf = LGBMClassifier()

        lgb_clf.fit(train_light, activity_train)

        y_pred = lgb_clf.predict(test_light)
        predictions = [round(value) for value in y_pred]

        lgb_accuracy = accuracy_score(activity_test, predictions)

        cf_lgb = confusion_matrix(activity_test, predictions)
        report_lgb = classification_report(activity_test, predictions)
        MCC_lgb = matthews_corrcoef(activity_test, predictions)

        return lgb_accuracy, cf_lgb, report_lgb, MCC_lgb


    @staticmethod
    def GBR_optimization(train, bindingEnergy_train, test, bindingEnergy_test):
        '''
        Parameter optimization for the Gradient Boosting Regressor.
        '''

        lr_list = [0.01, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
        R2_train = []
        R2_test = []
        best_test_R2 = None
        best_train_R2 = None
        best_lr = None
        for learning_rate in lr_list:
            gb_reg = GradientBoostingRegressor(n_estimators=10000, loss='ls', learning_rate=learning_rate, min_samples_split=3, max_depth=7, random_state=0)
            gb_reg.fit(train, bindingEnergy_train)

            R2_train.append(gb_reg.score(train, bindingEnergy_train))
            R2_test.append(gb_reg.score(test, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def GBR(best_lr):
        '''
        GradientBoostingRegressor algorithm.
        '''

        epoch_pred = []

        gb_reg = GradientBoostingRegressor(n_estimators=10000, loss='ls', learning_rate=best_lr, min_samples_split=3, max_depth=7, random_state=0)
        gb_reg.fit(train, bindingEnergy_train)


        for j, pred in enumerate(gb_reg.staged_predict(test)):
            epoch_pred.append((j,mean_squared_error(bindingEnergy_test, pred)))

        predictions = gb_reg.predict(test)
        R2_test = gb_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, epoch_pred


    @staticmethod
    def XGBR_optimization(train, bindingEnergy_train, test, bindingEnergy_test):
        '''
        Parameter optimization for the Xtreme Gradient Boosting Regressor.
        '''

        lr_list = [0.01, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
        R2_train = []
        R2_test = []
        best_test_R2 = None
        best_train_R2 = None
        best_lr = None
        for learning_rate in lr_list:
            xgb_reg = XGBRegressor(n_estimators=10000, learning_rate=learning_rate, random_state=0)
            xgb_reg.fit(train, bindingEnergy_train)

            R2_train.append(xgb_reg.score(train, bindingEnergy_train))
            R2_test.append(xgb_reg.score(test, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def XGBR(best_lr):
        '''
        XtremeGradientBoostingRegressor algorithm.
        '''

        epoch_pred = []

        xgb_reg =  XGBRegressor(n_estimators=10000, learning_rate=best_lr, random_state=0)
        xgb_reg.fit(train, bindingEnergy_train)


        predictions = xgb_reg.predict(test)
        R2_test = xgb_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE

    @staticmethod
    def LGBR_optimization(train, bindingEnergy_train, test, bindingEnergy_test):
        '''
        Parameter optimization for the Light Gradient Boosting Regressor.
        '''

        lr_list = [0.01, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
        R2_train = []
        R2_test = []
        best_test_R2 = None
        best_train_R2 = None
        best_lr = None
        for learning_rate in lr_list:
            lgb_reg = LGBMRegressor(n_estimators=10000, learning_rate=learning_rate, max_depth=7, random_state=0)
            lgb_reg.fit(train, bindingEnergy_train)

            R2_train.append(lgb_reg.score(train, bindingEnergy_train))
            R2_test.append(lgb_reg.score(test, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def LGBR(best_lr):
        '''
        LightGradientBoostingRegressor algorithm.
        '''

        lgb_reg = LGBMRegressor(n_estimators=10000, learning_rate=best_lr, max_depth=7, random_state=0)
        lgb_reg.fit(train, bindingEnergy_train)

        predictions = lgb_reg.predict(test)
        R2_test = lgb_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE

    @staticmethod
    def MLPR(train, bindingEnergy_train, test, bindingEnergy_test):
        '''
        Multi-layer Perceptor Regressor algorithm.
        '''

        parameter_space = {
                            'hidden_layer_sizes': [(100,)],
                            'activation': ['relu'],
                            'solver': ['adam'],
                            'alpha': [0.25],
                            'learning_rate': ['constant'],
}


        mlpr = MLPRegressor(random_state=1, max_iter=500)
        mlpr_reg = GridSearchCV(mlpr, parameter_space, cv=5)
        mlpr_reg.fit(train, bindingEnergy_train)

        best_params = mlpr_reg.best_params_

        predictions = mlpr_reg.predict(test)
        R2_test = mlpr_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, best_params



def main():

    start = time.time()

    path, test_size, seed, classification, regression, data_scaling = parseArg()

    LM = LearningModels()
    S = Scaler()

    global train
    global test
    global activity_train
    global activity_test
    global bindingEnergy_train
    global bindingEnergy_test

    data = read_graph_data(path=path)
    data = dataset_preparation(data=data)

    if classification:
        train, activity_train, test, activity_test = data_splitting_classification(data=data, test_size=test_size, seed=seed)
        best_lr, best_train_acc, best_test_acc = LM.GBC_optimization(train=train, activity_train=activity_train, test=test, activity_test=activity_test)
        cf, report, MCC = LM.GBC(best_lr=best_lr)
        XGB_accuracy, cf_XGB, report_XGB, MCC_XGB = LM.XGBoost(train=train, activity_train=activity_train, test=test, activity_test=activity_test)
        lgb_accuracy, cf_lgb, report_lgb, MCC_lgb = LM.LigthGB(train=train, activity_train=activity_train, test=test, activity_test=activity_test)

        print("#####################################")
        print("Gradient Boosting Classifier\n")
        print("Best learning rate: {}".format(best_lr))
        print("Best train accuracy: {}".format(best_train_acc))
        print("Best test accuracy: {}\n".format(best_test_acc))
        print("Confusion matrix")
        print(cf)
        print("\nClassification report:")
        print(report)
        print("\nMCC: {}".format(MCC))
        print("\n#####################################")
        print("XGBoost Classifier\n")
        print("XGB accuracy: {}\n".format(XGB_accuracy))
        print("Confusion matrix")
        print(cf_XGB)
        print("\nClassification report:")
        print(report_XGB)
        print("\nMCC: {}".format(MCC_XGB))
        print("\n#####################################")
        print("LightGBoost Classifier\n")
        print("LGBMC accuracy: {}\n".format(lgb_accuracy))
        print("Confusion matrix")
        print(cf_lgb)
        print("\nClassification report:")
        print(report_lgb)
        print("\nMCC: {}".format(MCC_lgb))

    elif regression:

        if data_scaling:

            scaling_functions = [(S.standard_scaler, "standard"), (S.min_max_scaler, 'min max'), (S.power_scaler, 'power'), (S.quantile_scaler_uniform, 'uniform'), (S.quantile_scaler_gauss, 'quantile Gauss'), (S.robust_scaler, 'robust'), (S.max_abs_scaler, 'max abs')]

            for scale in scaling_functions:

                train, bindingEnergy_train, test, bindingEnergy_test, ligands_test = data_splitting_regression(data=data, test_size=test_size, seed=seed)

                print("---------------------------------------")
                print("\nApplying {} scaler\n".format(scale[1]))

                try:
                    train, test = scale[0](train=train, test=test)

                except ValueError:
                    print("{} scaler cannot be applied to this data values.\n".format(scale[1]))
                    continue

                best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
                pred, R2_test, MSE = LM.GBR(best_lr=best_lr)
                MLPR_pred, MLPR_R2_test, MLPR_MSE, params = LM.MLPR(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)

                print("#####################################")
                print("Gradient Boosting Regressor\n")
                print("Best learning rate: {}".format(best_lr))
                print("Best train R2: {}".format(best_train_R2))
                print("Best test R2: {}\n".format(best_test_R2))
                print("Mean squared error: {}\n".format(MSE))
                print("#####################################")
                print("Multi-layer Perceptron Regressor\n")
                print("Best learning rate: {}".format(MLPR_best_lr))
                print("Best train R2: {}".format(MLPR_best_train_R2))
                print("Best parameters: {}".format(params))
                print("Best test R2: {}\n".format(MLPR_R2_test))
                print("Mean squared error: {}\n".format(MLPR_MSE))

        else:
            train, bindingEnergy_train, test, bindingEnergy_test, ligands_test = data_splitting_regression(data=data, test_size=test_size, seed=seed)
            best_lr, best_train_R2, best_test_R2 = LM.GBR_optimization(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
            pred, R2_test, MSE, epoch_pred = LM.GBR(best_lr=best_lr)
            XGBR_best_lr, XGBR_best_train_R2, XGBR_best_test_R2 = LM.XGBR_optimization(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
            XGBR_pred, XGBR_R2_test, XGBR_MSE = LM.XGBR(best_lr=XGBR_best_lr)
            LGBR_best_lr, LGBR_best_train_R2, LGBR_best_test_R2 = LM.LGBR_optimization(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)
            LGBR_pred, LGBR_R2_test, LGBR_MSE = LM.LGBR(best_lr=LGBR_best_lr)
            MLPR_pred, MLPR_R2_test, MLPR_MSE, params = LM.MLPR(train=train, bindingEnergy_train=bindingEnergy_train, test=test, bindingEnergy_test=bindingEnergy_test)

            print("#####################################")
            print("Gradient Boosting Regressor\n")
            print("Best learning rate: {}".format(best_lr))
            print("Best train R2: {}".format(best_train_R2))
            print("Best test R2: {}\n".format(best_test_R2))
            print("Mean squared error: {}\n".format(MSE))
            print("Predictions for epoch: \n")
            print("Delta MSE: {}\n".format(epoch_pred[-1][1] - epoch_pred[0][1]))
            for i in range(0, len(epoch_pred), 500):
                print("MSE at epoch {}: {}".format(epoch_pred[i][0], epoch_pred[i][1]))

            print("#####################################")
            print("Xtreme Gradient Boosting Regressor\n")
            print("Best learning rate: {}".format(XGBR_best_lr))
            print("Best train R2: {}".format(XGBR_best_train_R2))
            print("Best test R2: {}\n".format(XGBR_best_test_R2))
            print("Mean squared error: {}\n".format(XGBR_MSE))
            print("Predictions for epoch: \n")
            print("Delta MSE: {}\n".format(epoch_pred[-1][1] - epoch_pred[0][1]))
            for i in range(0, len(epoch_pred), 500):
                print("MSE at epoch {}: {}".format(epoch_pred[i][0], epoch_pred[i][1]))


            print("#####################################")
            print("Light Gradient Boosting Regressor\n")
            print("Best learning rate: {}".format(LGBR_best_lr))
            print("Best train R2: {}".format(LGBR_best_train_R2))
            print("Best test R2: {}\n".format(LGBR_best_test_R2))
            print("Mean squared error: {}\n".format(LGBR_MSE))

            print("#####################################")
            print("Multi-layer Perceptron Regressor\n")
            print("Best parameters: {}".format(params))
            print("Best test R2: {}\n".format(MLPR_R2_test))
            print("Mean squared error: {}\n".format(MLPR_MSE))
            print("MPL predicctions: {}\n".format(MLPR_pred))

    end = time.time()

    print("\nRun time: {} seconds.".format(end - start))


if __name__ == "__main__":
    main()
