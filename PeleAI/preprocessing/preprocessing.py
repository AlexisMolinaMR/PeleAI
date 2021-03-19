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
