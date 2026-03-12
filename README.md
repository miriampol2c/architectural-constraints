# Adaptive and Spandrel-like Constraints at Functional Sites in Protein Folds
📄 Read the full [preprint](https://www.biorxiv.org/content/10.64898/2026.02.09.704872v1) on bioRxiv.

Our methodology (c) combines reverse-folding and structure-prediction methods (a) with sequence analysis and local energetic frustration (b) to address several long-standing challenges in molecular biology:
- How protein sequence diversity maps onto a conserved three-dimensional fold.
- Which sequence elements primarily support structural integrity versus molecular function.
- Whether protein folds are shaped predominantly by evolutionary selection or by fundamental physical constraints.
- How local energetic frustration persists within globally minimally frustrated protein energy landscapes.
<img width="12503" height="11189" alt="final" src="https://github.com/user-attachments/assets/5bf8c4a5-03ec-4e05-a320-240a569a1232" />

### Main data
All input/output data needed to reproduce the main results of this article as well as the intermediate analyses are available at this ZENODO repository [TBC].

### Retrained ProteinMPNN versions
Customised datasets and models weights are available at this ZENODO repository [TBC].
|                        | PDB Hydrolases | PDB AllEnzymes | CLEAN Hydrolases | CLEAN AllEnzymes |
|:----------------------:|:--------------:|:--------------:|:----------------:|:----------------:|
| Training num epochs    |      120       |      120       |       120        |       120        |
| Backbone noise training|      0.2       |      0.2       |       0.2        |       0.2        |
| Training top_k         |       48       |       48       |        48        |        48        |
| Train dataset clusters |     22517      |     19113      |      21908       |      17362       |
| Valid datasets clusters|      1349      |      1167      |       1311       |       1042       |
| Test datasets clusters |      1424      |      1212      |       1378       |       1070       |
| Best model @ epoch    |       111      |       116      |        112       |        118       |

### Code
This GitHub repository contains the code for:
- #### Retrain and evaluate ProteinMPNN versions 
([train-and-eval-code](https://github.com/miriampol2c/architectural-constraints/tree/main/train-and-eval-code))

Needed: Python>=3.0, PyTorch, Numpy.

The multi-chain training data (16.5 GB, PDB biounits, 2021 August 2) can be downloaded with:
```
wget https://files.ipd.uw.edu/pub/training_sets/pdb_2021aug02.tar.gz
```
To extract:
```
tar -xvzf pdb_2021aug02.tar.gz
```
To retrain the model:
```
python3 train.py --data_path="path/of/extracted/dataset" --list_path="path/of/customised/data.csv"
```
For inference and evaluation:
```
python3 eval.py --data_path="path/of/extracted/dataset" --list_path="path/of/customised/data.csv" --model_path="path/to/your/new/weights.pt"
```
After Execution of each command read the results in the training.log or testing.log files in the 'loggs/' directory.
- #### Preliminary design pipelines applicable to any protein family
([generation-pipelines](https://github.com/miriampol2c/architectural-constraints/tree/main/main-code))

Currently both under development for personal usage, not yet optimized for end-user usage. Just possible to reuse individual scripts and functions.
- #### Code for the main figures included in the manuscript
