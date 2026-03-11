# Adaptive and Spandrel-like Constraints at Functional Sites in Protein Folds
📄 Read the full [preprint](https://www.biorxiv.org/content/10.64898/2026.02.09.704872v1) on bioRxiv.

Here we combined reverse-folding and structure-prediction methods with sequence analysis and local energetic frustration to address several long-standing challenges in molecular biology:
- How protein sequence diversity maps onto a conserved three-dimensional fold.
- Which sequence elements primarily support structural integrity versus molecular function.
- Whether protein folds are shaped predominantly by evolutionary selection or by fundamental physical constraints.
- How local energetic frustration persists within globally minimally frustrated protein energy landscapes.
<img width="12503" height="11189" alt="final" src="https://github.com/user-attachments/assets/5bf8c4a5-03ec-4e05-a320-240a569a1232" />

### Retrained ProteinMPNN versions
Models weights are available at this ZENODO repository [TBC].
|                        | PDB Hydrolases | PDB AllEnzymes | CLEAN Hydrolases | CLEAN AllEnzymes |
|:----------------------:|:--------------:|:--------------:|:----------------:|:----------------:|
| Training num epochs    |      120       |      120       |       120        |       120        |
| Backbone noise training|      0.2       |      0.2       |       0.2        |       0.2        |
| Training top_k         |       48       |       48       |        48        |        48        |
| Train dataset clusters |     22517      |     19113      |      21908       |      17362       |
| Valid datasets clusters|      1349      |      1167      |       1311       |       1042       |
| Test datasets clusters |      1424      |      1212      |       1378       |       1070       |
| Best models @ epoch    |       111      |       116      |        112       |        118       |

### Main data
All input/output data needed to reproduce the main results of this article as well as the intermediate analysis are available at this ZENODO repository [TBC].

### Code
This GitHub repository contains:
- Code to retrain and evaluate ProteinMPNN versions.
- Preliminar code to apply our design pipeline on different protein families (both on development).
- Code for the main figures included in the manuscript. 
