import pandas as pd
import glob
import matplotlib.pyplot as plt # to seaborn
import os
import sys

# to argparse

path = sys.argv[1] # path by sys argv for now
metric = sys.argv[2] # plot metric vs B.E.
show = sys.argv[3] # show if set to 'show', save if set to 'save'
save_path = sys.argv[4] # path to save the plots

fields = ['hbond_H_val_690', 'RMSD_ligand', 'docking_com_dist', 'Binding Energy', 'trajectory'] # to be given as input (argparse)

for filename in glob.glob(os.path.join(path, '*.csv')):
    pele_out_summary = pd.read_csv(filename, sep=';', usecols = fields)
    pele_out_summary.rename({'Binding Energy': 'Binding_Energy'}, axis = 'columns', inplace = True)

    if metric == "hbond":
        plt.scatter(pele_out_summary.hbond_H_val_690, pele_out_summary.Binding_Energy)
        plt.xlabel("H-bond")
    elif metric == "rmsd":
        plt.scatter(pele_out_summary.RMSD_ligand, pele_out_summary.Binding_Energy)
        plt.xlabel("RMSD")
    elif metric == "docking":
        plt.scatter(pele_out_summary.docking_com_dist, pele_out_summary.Binding_Energy)
        plt.xlabel("Docking distance")

    plt.ylabel("Binding Energy")
    plt.title(pele_out_summary.trajectory[0].split('/')[0])

    if show == 'show':
        plt.show()
    else:
        plt.savefig(save_path + '/' + pele_out_summary.trajectory[0].split('/')[0] + '_' + metric + "_pele_plot.png")

    plt.close()
