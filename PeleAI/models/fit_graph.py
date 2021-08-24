import re
import os
import warnings

import tensorflow as tf

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, matthews_corrcoef, mean_squared_error, mean_absolute_error
from sklearn.model_selection import StratifiedKFold, train_test_split, GridSearchCV, RandomizedSearchCV
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn import feature_selection
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier, XGBRegressor
from lightgbm import LGBMClassifier, LGBMRegressor

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras import optimizers
from keras.layers import Dropout
from tensorflow.keras.optimizers import Nadam

from utils.utils import r_squared

os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
warnings.filterwarnings("ignore")


class LearningModels():
    '''
    This class contains all learning algorithms used for fitting the data.
    '''

    global train

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
            gb_clf = GradientBoostingClassifier(
                n_estimators=20, learning_rate=learning_rate, max_features=2, max_depth=2, random_state=0)
            gb_clf.fit(train, activity_train)

            accuracies_train.append(gb_clf.score(train, activity_train))
            accuracies_test.append(gb_clf.score(test, activity_test))

        best_test_acc = max(accuracies_test)
        best_train_acc = accuracies_train[accuracies_test.index(best_test_acc)]
        best_lr = lr_list[accuracies_test.index(best_test_acc)]

        return best_lr, best_train_acc, best_test_acc

    @staticmethod
    def GBC(train, activity_train, test, activity_test, best_lr):
        '''
        GradientBoostingClassifier algorithm.
        '''

        gb_clf2 = GradientBoostingClassifier(
            n_estimators=20, learning_rate=best_lr, max_features=2, max_depth=2, random_state=0)
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

        train_light = train.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
        test_light = test.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))

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
            gb_reg = GradientBoostingRegressor(
                n_estimators=10000, loss='ls', learning_rate=learning_rate, min_samples_split=3, max_depth=7, random_state=0)
            gb_reg.fit(train, bindingEnergy_train)

            R2_train.append(gb_reg.score(train, bindingEnergy_train))
            R2_test.append(gb_reg.score(test, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def GBR(train, bindingEnergy_train, test, bindingEnergy_test, best_lr):
        '''
        GradientBoostingRegressor algorithm.
        '''

        epoch_pred = []

        gb_reg = GradientBoostingRegressor(
            n_estimators=10000, loss='ls', learning_rate=best_lr, min_samples_split=3, max_depth=7, random_state=0)
        gb_reg.fit(train, bindingEnergy_train)

        for j, pred in enumerate(gb_reg.staged_predict(test)):
            epoch_pred.append((j, mean_squared_error(bindingEnergy_test, pred)))

        predictions = gb_reg.predict(test)
        R2_test = gb_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)
        MAE = mean_absolute_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, MAE, epoch_pred

    @staticmethod
    def XGBR_optimization(train, bindingEnergy_train, test, bindingEnergy_test, cpus):
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
            xgb_reg = XGBRegressor(
                n_estimators=10000, learning_rate=learning_rate, random_state=0, n_jobs=cpus)
            xgb_reg.fit(train, bindingEnergy_train)

            R2_train.append(xgb_reg.score(train, bindingEnergy_train))
            R2_test.append(xgb_reg.score(test, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def XGBR(train, bindingEnergy_train, test, bindingEnergy_test, best_lr, cpus):
        '''
        XtremeGradientBoostingRegressor algorithm.
        '''

        epoch_pred = []

        xgb_reg = XGBRegressor(n_estimators=10000, learning_rate=best_lr,
                               random_state=0, n_jobs=cpus)
        xgb_reg.fit(train, bindingEnergy_train)

        predictions = xgb_reg.predict(test)
        R2_test = xgb_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)
        MAE = mean_absolute_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, MAE

    @staticmethod
    def LGBR_optimization(train, bindingEnergy_train, test, bindingEnergy_test, cpus):
        '''
        Parameter optimization for the Light Gradient Boosting Regressor.
        '''

        lr_list = [0.01, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
        R2_train = []
        R2_test = []
        best_test_R2 = None
        best_train_R2 = None
        best_lr = None

        train_light = train.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
        test_light = test.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))

        for learning_rate in lr_list:
            lgb_reg = LGBMRegressor(
                n_estimators=10000, learning_rate=learning_rate, max_depth=15, random_state=0, n_jobs=cpus)
            lgb_reg.fit(train_light, bindingEnergy_train)

            R2_train.append(lgb_reg.score(train_light, bindingEnergy_train))
            R2_test.append(lgb_reg.score(test_light, bindingEnergy_test))

        best_test_R2 = max(R2_test)
        best_train_R2 = R2_train[R2_test.index(best_test_R2)]
        best_lr = lr_list[R2_test.index(best_test_R2)]

        return best_lr, best_train_R2, best_test_R2

    @staticmethod
    def LGBR(train, bindingEnergy_train, test, bindingEnergy_test, best_lr, cpus):
        '''
        LightGradientBoostingRegressor algorithm.
        '''

        train_light = train.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
        test_light = test.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))

        lgb_reg = LGBMRegressor(n_estimators=10000, learning_rate=best_lr,
                                max_depth=15, random_state=0, n_jobs=cpus)

        lgb_reg.fit(train_light, bindingEnergy_train)

        predictions = lgb_reg.predict(test_light)
        R2_test = lgb_reg.score(test_light, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)
        MAE = mean_absolute_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, MAE

    @staticmethod
    def MLPR(train, bindingEnergy_train, test, bindingEnergy_test, cpus):
        '''
        Multi-layer Perceptron Regressor algorithm.
        '''

        parameter_space = {
            'hidden_layer_sizes': [(100,)],
            'activation': ['relu'],
            'solver': ['adam'],
            'alpha': [0.25],
            'learning_rate': ['constant'],
        }

        mlpr = MLPRegressor(random_state=1, max_iter=500)
        mlpr_reg = GridSearchCV(mlpr, parameter_space, cv=5, n_jobs=cpus)
        mlpr_reg.fit(train, bindingEnergy_train)

        best_params = mlpr_reg.best_params_

        predictions = mlpr_reg.predict(test)
        R2_test = mlpr_reg.score(test, bindingEnergy_test)
        MSE = mean_squared_error(bindingEnergy_test, predictions)
        MAE = mean_absolute_error(bindingEnergy_test, predictions)

        return predictions, R2_test, MSE, MAE, best_params

    @staticmethod
    def FFNN_profile(train, val, test, learning_rate):
        '''
        '''

        NN_model = Sequential()

        NN_model.add(Dense(2048, kernel_initializer='normal',input_dim = train.shape[1], activation=tf.keras.layers.LeakyReLU()))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))

        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))

        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2048,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(1024,activation=tf.keras.layers.LeakyReLU(), use_bias=True))

        NN_model.add(Dense(1, activation='linear'))

        opt = tf.keras.optimizers.Nadam(learning_rate=learning_rate)

        NN_model.compile(loss=tf.keras.losses.Huber(), optimizer=opt, metrics=['mse','mae',r_squared])
        NN_model.summary()

        return NN_model

    @staticmethod
    def FFNN_clustering(train, val, test, learning_rate):
        '''
        '''

        NN_model = Sequential()

        NN_model.add(Dense(64, kernel_initializer='normal',input_dim = train.shape[1], activation=tf.keras.layers.LeakyReLU()))
        NN_model.add(Dense(64,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(32,activation=tf.keras.layers.LeakyReLU(), use_bias=True))

        NN_model.add(Dense(16,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(8,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(4,activation=tf.keras.layers.LeakyReLU(), use_bias=True))
        NN_model.add(Dense(2,activation=tf.keras.layers.LeakyReLU(), use_bias=True))

        NN_model.add(Dense(1, activation='linear'))

        opt = tf.keras.optimizers.Nadam(learning_rate=learning_rate)

        NN_model.compile(loss=tf.keras.losses.Huber(), optimizer=opt, metrics=['mse','mae',r_squared])
        NN_model.summary()

        return NN_model
