import pandas as pd
import glob
import matplotlib.pyplot as plt
import os
import sys


path = sys.argv[1] # path by sys argv for now
metric = sys.argv[2] # plot metric vs B.E.
show = sys.argv[3] # show if set to 'show', save if set to 'save'

fields = ['hbond_H_val_690', 'RMSD_ligand', 'docking_com_dist', 'Binding Energy'] # to be given as input

for filename in glob.glob(os.path.join(path, 'summary.csv')):
    pele_out_summary = pd.read_csv(filename, sep=';', usecols = fields)
    pele_out_summary.rename({'Binding Energy': 'Binding_Energy'}, axis = 'columns', inplace = True)

    foldername = os.path.basename(path)

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
    plt.title(foldername)

    if show == 'show':
        plt.show()
    else:
        plt.savefig(foldername + '_' + metric + "_pele_plot.png")
