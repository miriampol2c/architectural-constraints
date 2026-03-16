# Adaptive and Spandrel-like Constraints at Functional Sites in Protein Folds

📄 Read the full [preprint](https://www.biorxiv.org/content/10.64898/2026.02.09.704872v1) on bioRxiv [![DOI](https://img.shields.io/badge/DOI-10.1101%2F2026.02.09.704872-blue)](https://www.biorxiv.org/content/10.64898/2026.02.09.704872v1)

## Project Overview
Understanding the relationships among amino acid sequences, structures and functions in proteins and how they evolve, remains a central challenge in molecular biology. It is still unclear which sequence elements differentially contribute to structural integrity or molecular function. Even more, there are ongoing debates on whether protein folds emerge as a result of evolution or as a consequence of physical laws. The energy landscapes theory states that proteins are minimally frustrated systems, i.e. they fold by minimising their energetic conflicts. However, some local frustration, believed to be selected for functional reasons, remains in the native state of proteins. Here, we combine reverse folding and structure prediction methods with sequence and local frustration analysis to address the aforementioned ideas. We found that reverse folding techniques are unable to erase evolutionary conserved frustration from certain residues, even when detrimental for structural integrity. We propose that certain frustration hotspots behave like architectural spandrels, not directly shaped by selection but emerging from physical constraints in protein folds which evolution can later co-opt for function. Our results provide a new perspective revealing how sequence variation and functional specificity could evolve from evolutionary, structural and biophysical constraints.
<img width="12503" height="11189" alt="final" src="https://github.com/user-attachments/assets/5bf8c4a5-03ec-4e05-a320-240a569a1232" />
Figure: Exploration of frustration patterns and biophysical constraints of reverse-folded sequences. (a) Reverse-folded sequences are first generated and structurally predicted to confirm that they match a given backbone. (b) Local frustration is then computed for each design using FrustratometeR (Rausch et al., 2021) and integrated across all designs with FrustraEvo (Parra et al., 2024), that measures frustration conservation in the columns of the MSA. Because residue identities may vary among designs, the frustration state at each position may also differ across models. By combining frustration information from all sequences and their MSA, FrustraEvo highlights conserved and variable frustration patterns, revealing which positions share similar energetic behavior and which can tolerate alternative residues. Considering the positions remodelled versus those that remain largely invariant, can help to identify regions constrained by folding, stability, or function. (c) Overview of the proposed workflow for analysing frustration and biophysical constraints in any target or protein family. Combining (a) and (b) we implement two complementary strategies. The Single Target approach designs multiple sequences for one target, their predicted models are then analysed collectively to assess frustration conservation across the entire design set. In contrast, the Family Target approach evaluates a reverse-folding algorithm on a set of homologous structures. For each family member, a best-scoring design (e.g., highest pLDDT) is selected, building a synthetic designed family. FrustraEvo then quantifies frustration conservation across this designed family, providing deeper coverage of sequence space and offering insight into how reverse-folding algorithms handle conserved structural and functional constraints.

## Repository Structure
This GitHub repository contains the code we used for:
### (1) Run CLEAN function prediction
First be sure you have installed CLEAN as documented originally in [here](https://github.com/tttianhao/CLEAN/tree/main).
Then follow the instructions or check our scripts directly [CLEAN_scripts](https://github.com/miriampol2c/architectural-constraints/tree/main/CLEAN_scripts).
To generate embeddings:
```
python3 CLEAN_compute_emb.py --fasta_data your_fasta_name
```
For inference:
```
python3 CLEAN_infer_fasta.py --fasta_data your_fasta_name
```
After execution check the ‘/CLEAN/Aapp/results/’ directory, as the CLEAN software is designed to save the predictions results there as a csv file. 

### (2) Retrain and evaluate ProteinMPNN versions 
See [train-and-eval-code](https://github.com/miriampol2c/architectural-constraints/tree/main/train-and-eval-code)

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
After execution of each command read the results in the training.log or testing.log files in the 'loggs/' directory.
  
### (3) Run design pipelines applicable to any protein family
See [generation-pipelines](https://github.com/miriampol2c/architectural-constraints/tree/main/main-code)

Currently both under development for personal usage, not yet optimized for end-user usage. Just possible to reuse individual scripts and functions.

### (4) Code for the main figures included in the manuscript
See [code-main-figures](https://github.com/miriampol2c/architectural-constraints/tree/main/code-main-figures)

We provided the code to reproduce Multiple Sequence Frustration Alignment (MSFA) and Consensus MSFA visualizations (considering the frustration conservation threshold for coloring instead of the median index per position), following latest frustration analyses Freiberger et al., 2023](https://www.nature.com/articles/s41467-023-43801-2); as well as the code to generate contacts networks by using frustration results.

## Data Availability
### Results
All input/output data needed to reproduce the main results of this article as well as the intermediate analyses are available at this ZENODO repository [![](https://img.shields.io/badge/DOI-10.1038/zenodo.18922047-blue)](https://doi.org/10.5281/zenodo.18922047)

The results data structure is generally described there, check the [repo](https://doi.org/10.5281/zenodo.18922047).
### Retrained ProteinMPNN versions
Customised datasets and models weights are available at this ZENODO repository [![](https://img.shields.io/badge/DOI-10.1038/zenodo.18957473-blue)](https://doi.org/10.5281/zenodo.18957473)
|                        | PDB Hydrolases | PDB AllEnzymes | CLEAN Hydrolases | CLEAN AllEnzymes |
|:----------------------:|:--------------:|:--------------:|:----------------:|:----------------:|
| Training num epochs    |      120       |      120       |       120        |       120        |
| Backbone noise training|      0.2       |      0.2       |       0.2        |       0.2        |
| Training top_k         |       48       |       48       |        48        |        48        |
| Train dataset clusters |     22517      |     19113      |      21908       |      17362       |
| Valid datasets clusters|      1349      |      1167      |       1311       |       1042       |
| Test datasets clusters |      1424      |      1212      |       1378       |       1070       |
| Best model @ epoch    |       111      |       116      |        112       |        118       |

## Citation
### Manuscript
  [Citation information to be added upon publication] - for now:
   ```
  Poley-Gil, M., Fernandez-Martin, M., Banka, A., Heinzinger, M., Rost, B., Valencia, A., & Parra, R. G. (2026). Adaptive and Spandrel-like Constraints at Functional Sites in Protein Folds. bioRxiv, 2026-02.
   ```
### Data
  ```
  Poley-Gil, M. (2026). Architectural-constraints-main-data [Data set]. Zenodo. https://doi.org/10.5281/zenodo.18922047
  ```
### Models
  ```
  Poley-Gil, M. (2026). Retrained-ProteinMPNN-versions. Zenodo. https://doi.org/10.5281/zenodo.18957473
  ```
