import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler, RobustScaler, Normalizer, QuantileTransformer, PowerTransformer


def dataset_preparation(data):
    '''
    Prepare dataframe for further data preparation.
    '''

    data.reset_index(drop=True, inplace=True)
    data.fillna(0, inplace=True)

    return data


def data_splitting_classification(data, test_size, seed):
    '''
    Split the date in train and test sets if classification problem.
    '''

    train, test = train_test_split(data, test_size=test_size, random_state=seed)

    activity_train = train['activity']
    activity_test = test['activity']

    train.drop(['ligand', 'activity'], axis=1, inplace=True)
    test.drop(['ligand', 'activity'], axis=1, inplace=True)

    return train, activity_train, test, activity_test


def data_splitting_regression(data, test_size, seed):
    '''
    Split the date in train and test sets if regression problem.
    '''

    train, test = train_test_split(data, test_size=test_size, random_state=seed)

    bindingEnergy_train = train['bindingEnergy']
    bindingEnergy_test = test['bindingEnergy']

    ligands_test = test['ligand']

    train.drop(['ligand', 'bindingEnergy'], axis=1, inplace=True)
    test.drop(['ligand', 'bindingEnergy'], axis=1, inplace=True)

    return train, bindingEnergy_train, test, bindingEnergy_test, ligand_test

def data_splitting_ffnn(data, seed):
    '''
    Split the date in train and test sets if regression problem.
    '''

    train, val, test = np.split(data.sample(frac=1, random_state=seed), [int(.6*len(data)), int(.8*len(data))])

    if 'bindingEnergy' in list(data.columns):

        bindingEnergy_train = train['bindingEnergy']
        bindingEnergy_test = test['bindingEnergy']
        bindingEnergy_val = val['bindingEnergy']

        train.drop(['bindingEnergy'],axis=1, inplace=True)
        test.drop(['bindingEnergy'], axis=1, inplace=True)
        val.drop(['bindingEnergy'], axis=1, inplace=True)


    else:

        bindingEnergy_train = train['IC50']
        bindingEnergy_test = test['IC50']
        bindingEnergy_val = val['IC50']

        train.drop(['IC50'],axis=1, inplace=True)
        test.drop(['IC50'], axis=1, inplace=True)
        val.drop(['IC50'], axis=1, inplace=True)


    if 'ligandRMSD' in data.columns:

        ligandRMSD_test = test['ligandRMSD']

        train.drop(['ligandRMSD'],axis=1, inplace=True)
        test.drop(['ligandRMSD'], axis=1, inplace=True)
        val.drop(['ligandRMSD'], axis=1, inplace=True)

    else:

        ligandRMSD_test = None

    ligand_test = test['ligand']

    train.drop(['ligand'],axis=1, inplace=True)
    test.drop(['ligand'], axis=1, inplace=True)
    val.drop(['ligand'], axis=1, inplace=True)

    return train, bindingEnergy_train, val, bindingEnergy_val, test, bindingEnergy_test, ligand_test, ligandRMSD_test


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
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def min_max_scaler(train, test):
        '''
        Transform features by scaling each feature to a given range. Set on default (0,1).
        '''

        scaler = MinMaxScaler()

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def min_max_scaler_ffnn(train, val, test):
        '''
        Transform features by scaling each feature to a given range. Set on default (0,1).
        '''

        scaler = MinMaxScaler()

        x = train.values
        x_scaled = scaler.fit_transform(x)
        train = pd.DataFrame(x_scaled)
        val_s = scaler.transform(val)
        val = val_s
        test_s = scaler.transform(test)
        test = test_s

        return train, val, test

    @staticmethod
    def max_abs_scaler(train, test):
        '''
        Scale each feature by its maximum absolute value.
        '''

        scaler = MaxAbsScaler()

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def robust_scaler(train, test):
        '''
        Scaler removes the median and scales the data according
        to the quantile range (defaults to IQR: Interquartile Range).
        '''

        scaler = RobustScaler()

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def normal_scaler(train, test):
        '''
        Normalize samples individually to unit norm.
        '''

        scaler = Normalizer()

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def quantile_scaler_gauss(train, test):
        '''
        Transforms the features to follow a normal distribution.
        '''

        scaler = QuantileTransformer(output_distribution='normal')

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def quantile_scaler_uniform(train, test):
        '''
        Transforms the features to follow a uniform distribution.
        '''

        scaler = QuantileTransformer(output_distribution='uniform')

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test

    @staticmethod
    def power_scaler(train, test):
        '''
        Apply a power transform featurewise to make data more Gaussian-like.
        '''

        scaler = PowerTransformer(method='yeo-johnson')

        train = scaler.fit_transform(train)
        test = scaler.transform(test)

        return train, test
