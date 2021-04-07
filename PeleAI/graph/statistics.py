import numpy as np


def laplacianStats(matrices, pose):
    '''
    Compute statistics out of Laplacian matrices
    '''

    laplacian_statistics = []

    print("\nComputing statistics for laplacian matrices...")

    laplacian_statistics.append(pose.split('.')[0])

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


def adjacencyStats(matrices, pose):
    '''
    Compute statistics out of adjacency matrices.
    '''

    adjacency_statistics = []

    print("\nComputing statistics for adjacency matrices...")

    adjacency_statistics.append(pose.split('.')[0])

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
                maximum = np.real(maximum)
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
