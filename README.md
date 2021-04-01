## Table of contents
* [PeleAI](#PeleAI)
* [Requirements](##Requirements)
* [Usage](#Usage)

# PeleAI

The goal of PeleAI-3D is to create target-specific scoring functions from simulation data obtained from molecular simulations methodologies such as Monte-Carlo or docking.

To do so an graph-based topological description of the binding site containing the ligand is computed, which serves later for fitting a model that predicts either the activitiy or the binding energy of the given pose. 

## Requirements

Install the conda version of ```rdkit``` and activate the evironment before installing the required packages:

```
conda create -c conda-forge -n my-rdkit-env rdkit
conda activate my-rdkit-env
```
Then, you can install:

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

The methodology can be run in three different fashions: _graph_, _fitting_ or _pipeline_. 

The first one, _graph_, is dedicated to generate statistics out from a graph topology of the ligand-binding site complex without fitting the data in any model. The second one, _fitting_, will provide several options to use the graph based statistics previously generated to fit either a regression or a classification task.
The last one, _pipeline_, provides end-to-end results out of the given input poses. It executes _graph_, whose output is passed to _fitting_ which will report the model results.

### Generate graph statistics (_graph_)

**Parameters**

```
#control file for PeleAI3D - Graph statistics

path: 
output: 
run_name: 
ligand_name: 
selection_radius: 
center: 
decay_function:
nodes: 
``` 

### Fitting graph statistics (_fitting_)

**Parameters**

```
#control file for PeleAI3D - Fitting a model

path_graph: 
output: 
test_size:
seed: 
task: 
cpus: 
scaler: 
algorithm: 
``` 

### Pipelining (_pipeline_)

**Parameters**

```
#control file for PeleAI3D - Pipeline
path: 
output: 
run_name: 
ligand_name: 
selection_radius: 
center: 
decay_function:
nodes: 
# ------------------------------------
#pipe: True
#target: 
# ------------------------------------
test_size: 
seed: 
task: 
cpus: 
scaler: 
algorithm: 
``` 

### Execution

You may pass the _input.yaml_ file to the ```peleAI3d.py``` as follows:

```
python3 peleAI3d.py input.yaml
```

The results will be written to the folder indicated in the ```output``` argument of the input file.
