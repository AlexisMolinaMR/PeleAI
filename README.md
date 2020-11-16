## Table of contents
* [PeleAI](#PeleAI)
* [Requirements](##Requirements)
* [Usage](##Usage)
* [To developers](##To-developers)

# PeleAI

Create scoring functions from simulation.

## Requirements

```
numpy==1.19.0
pandas==1.1.4
prody==1.11
sklearn==0.23
xgboost==1.3.0
lightgbm==3.1.0
```

`pip install -r requirements.txt`

## Usage

### Generate graph statistics

**Arguments**

`python3  graph_generator.py [-h] -i POSE -l LIGAND -r RADIUS [-gc] [-ds] (-ex | -lo) [-lp]`


**Basic usage (recomended)**

`python3 generate_graph.py -i path/to/pose.pdb -l LIG -r RADIUS -gc -ex`

You will need a bash script to run through several poses.

The output will be two csv files, one for Laplacian matrices (_L_stats_ex.csv_) and another for Adjacency matrices (_A_stats_ex.csv_). Use _L_stats_ex.csv_.

### Fit models

**Arguments**

`python3 -i /path/to/L_stats_ex.csv -t TEST_SIZE -s SEED (-c | -r) [-d]`

#### Classification

`python3 -i /path/to/L_stats_ex.csv -t TEST_SIZE -s SEED -c`

Classification will be performed using three gradient boosting regressors: Gradient Boosting Regressor, eXtreme Gradient Boosting Regressor and Light Gradient Boosting Regresssor.

##### Data scaling

If desired, all sklearn scalers will be fit and evaluated.

`python3 -i /path/to/L_stats_ex.csv -t TEST_SIZE -s SEED -c -d`

#### Regression

Not tested. Avoid for now.

## To developers

### Protein-ligand branch

Commit to this branch if developing for protein-ligand interactions.

### Protein-protein branch

Commit to this brach if developing for protein-protein interactions.
