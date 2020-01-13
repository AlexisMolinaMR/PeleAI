import pandas as pd
import glob
import matplotlib.pyplot as plt
import os
import sys


path = sys.argv[1] #path by sys argv for now

fields = ['hbond_H_val_690', 'Binding Energy'] # to be given as input

for filename in glob.glob(os.path.join(path, 'summary.csv')):
    pele_out_summary = pd.read_csv(filename, sep=';', usecols = fields)
    pele_out_summary.rename({'Binding Energy': 'Binding_Energy'}, axis = 'columns', inplace = True)
    plt.scatter(pele_out_summary.hbond_H_val_690, pele_out_summary.Binding_Energy)
    plt.show()
