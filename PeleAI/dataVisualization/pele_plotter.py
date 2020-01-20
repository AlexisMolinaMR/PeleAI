import pandas as pd
import glob
import matplotlib.pyplot as plt # to seaborn
import os
import sys
import plotly.express as px
import plotly.graph_objects as go
import math
# to argparse

path = sys.argv[1] # path by sys argv for now
metric = sys.argv[2] # plot metric vs B.E.
show = sys.argv[3] # show if set to 'show', save if set to 'save'
save_path = sys.argv[4] # path to save the plots
mode = sys.argv[5] # normal/interactive
filter = sys.argv[6] # select best %, NA to avoid

# e.g. python pele_plotter.py /home/alexis/Desktop/summaries hbond save /home/alexis/Desktop/plots interactive NA

fields = ['hbond_H_val_690', 'RMSD_ligand', 'docking_com_dist', 'Binding Energy', 'trajectory', 'epoch', 'numberOfAcceptedPeleSteps'] # to be given as input (argparse)

for filename in glob.glob(os.path.join(path, '*.csv')):
    pele_out_summary = pd.read_csv(filename, sep=';', usecols = fields)
    pele_out_summary.rename({'Binding Energy': 'Binding_Energy'}, axis = 'columns', inplace = True)

    pele_out_summary.numberOfAcceptedPeleSteps = pele_out_summary.numberOfAcceptedPeleSteps.astype(str)
    pele_out_summary.epoch = pele_out_summary.epoch.astype(str)


    if mode == "normal" and filter == 'NA':
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
            print("Writing image " + pele_out_summary.trajectory[0].split('/')[0] + ".png" + " to: " + save_path)

        plt.close()

    if mode == "normal" and filter != 'NA':
        filter_number = math.floor(len(pele_out_summary.Binding_Energy) * float(filter))
        filtered_dataframe = pele_out_summary.sort_values(by = "Binding_Energy", inplace = False)
        filtered_dataframe = filtered_dataframe.nlargest(filter_number, ['Binding_Energy'])
        if metric == "hbond":
            plt.scatter(filtered_dataframe.hbond_H_val_690, filtered_dataframe.Binding_Energy)
            plt.xlabel("H-bond")
        elif metric == "rmsd":
            plt.scatter(filtered_dataframe.RMSD_ligand, filtered_dataframe.Binding_Energy)
            plt.xlabel("RMSD")
        elif metric == "docking":
            plt.scatter(filtered_dataframe.docking_com_dist, filtered_dataframe.Binding_Energy)
            plt.xlabel("Docking distance")

        plt.ylabel("Binding Energy")
        plt.title(pele_out_summary.trajectory[0].split('/')[0] + " best {} poses".format(filter_number))

        if show == 'show':
            plt.show()
        else:
            plt.savefig(save_path + '/' + pele_out_summary.trajectory[0].split('/')[0] + '_' + metric + "_filtered_pele_plot.png")
            print("Writing image " + pele_out_summary.trajectory[0].split('/')[0] + "best {} poses.png to: ".format(filter_number) + save_path)

        plt.close()


    elif mode == "interactive" and filter == 'NA':
        print(pele_out_summary.trajectory.str.split("/"))
        if metric == "hbond":
            fig = go.Figure(data=go.Scatter(x=pele_out_summary.hbond_H_val_690,
                            y=pele_out_summary.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=pele_out_summary.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0],
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="H-bond distance",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )

        elif metric == "rmsd":
            fig = go.Figure(data=go.Scatter(x=pele_out_summary.RMSD_ligand,
                            y=pele_out_summary.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=pele_out_summary.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0],
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="RMSD",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )
        elif metric == "docking":
            fig = go.Figure(data=go.Scatter(x=pele_out_summary.docking_com_dist,
                            y=pele_out_summary.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=pele_out_summary.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0],
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="Docking distance",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )
        if show == 'show':
            fig.show()
        else:
            fig.write_image(save_path + '/' + pele_out_summary.trajectory[0].split('/')[0] + '_' + metric + "_interactive_pele_plot.png")
            print("Writing image " + pele_out_summary.trajectory[0].split('/')[0] + ".png" + " to: " + save_path)

    elif mode == "interactive" and filter != 'NA':
        filter_number = math.floor(len(pele_out_summary.Binding_Energy) * float(filter))
        filtered_dataframe = pele_out_summary.sort_values(by = "Binding_Energy", inplace = False)
        filtered_dataframe = filtered_dataframe.nlargest(filter_number, ['Binding_Energy'])
        if metric == "hbond":
            fig = go.Figure(data=go.Scatter(x=filtered_dataframe.hbond_H_val_690,
                            y=filtered_dataframe.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=filtered_dataframe.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0]+ " best {} poses".format(filter_number),
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="H-bond distance",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )

        elif metric == "rmsd":
            fig = go.Figure(data=go.Scatter(x=filtered_dataframe.RMSD_ligand,
                            y=filtered_dataframe.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=filtered_dataframe.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0]+ " best {} poses".format(filter_number),
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="RMSD",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )
        elif metric == "docking":
            fig = go.Figure(data=go.Scatter(x=filtered_dataframe.docking_com_dist,
                            y=filtered_dataframe.Binding_Energy,
                            mode = 'markers',
                            text = pele_out_summary.trajectory +'<br>'+ 'Epoch: '+ pele_out_summary.epoch +'<br>'+ 'Step: ' + pele_out_summary.numberOfAcceptedPeleSteps,
                            marker=dict(
                                size=10,
                                color=filtered_dataframe.Binding_Energy, #set color equal to a variable
                                colorscale='Viridis', # one of plotly colorscales
                                showscale=True)))
            fig.update_layout(title=pele_out_summary.trajectory[0].split('/')[0]+ " best {} poses".format(filter_number),
                              annotations=[
                                    go.layout.Annotation(
                                        x=0.5,
                                        y=-0.06,
                                        showarrow=False,
                                        text="Docking distance",
                                        xref="paper",
                                        yref="paper"
                                    ),
                                    go.layout.Annotation(
                                        x=-0.07,
                                        y=0.5,
                                        showarrow=False,
                                        text="Binding Energy",
                                        textangle=-90,
                                        xref="paper",
                                        yref="paper"
                                    )
                                ]
                            )
        if show == 'show':
            fig.show()
        else:
            fig.write_image(save_path + '/' + pele_out_summary.trajectory[0].split('/')[0] + '_' + metric + "_filtered_interactive_pele_plot.png")
            print("Writing image " + pele_out_summary.trajectory[0].split('/')[0] + " best {} poses.png to: ".format(filter_number) + save_path)

print("DONE!")