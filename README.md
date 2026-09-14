# H3-K27M mutation alters the dynamics of human hematopoietic stem cells and delays erythroid differentiation
This repository accompanies the article "H3-K27M mutation alters the dynamics of human hematopoietic stem cells and delays erythroid differentiation". It includes the original blood cell count tables from the xenotransplantation mouse model (XLSX-files "Original count table - EXP1 & EXP2") and the necessary code for the mathematical analysis.

Authors: Mia Brunetti<sup>1,2</sup>; Hassan Dakik<sup>3</sup>; Fatemeh Beigmohammadi<sup>2</sup>; Kolja Eppert<sup>3,4</sup>; Morgan Craig<sup>1,2</sup><br>
<sup>1</sup>Département de Mathématiques et de Statistiques, Université de Montréal, Montréal, Canada<br>
<sup>2</sup>Sainte-Justine University Hospital Research Center, Montréal, Canada<br>
<sup>3</sup>Department of Pediatrics, McGill University, Montréal, Canada<br>
<sup>4</sup>Research Institute of the McGill University Health Centre, Canada<br>

## Repository structure
| Folder | Description |
|:-----|:------------|
| Structural identifiability | MATLAB M-files necessary to create the ODE models stored in MAT-files and read by STRIKE-GOLDD. |
| Data | CSV-files containing the blood cell counts and Julia script initializing the data. |
| Models | Julia scripts for fitting hematopoiesis ODE models. |
| Parameter estimation | Adaptive simulated annealing (ASA) functions and commands for parameter estimation. |
| Parsimony analysis | Functions and commands for parsimony analysis. |
| TM-RWFS | Trajectory-Matching Random-Walk Feasibility Sampling (TM-RWFS) commands for parameter distribution. |
| Statistical test | Mann-Whitney U test on parameter distribution commands. |
| Results | CSV-files containing the results for parameter estimation, 95% credible intervals, and accepted trajectories for the selected models model. |
| Figure visualization | Functions, template, and commands for reproducing result figures. |

## Requirements
- MATLAB R2025b or later for the structural identifiability analysis.
- Julia 1.11.9 or later.
- For the Julia environment, the repository does not include a `Project.toml`. Install the required packages as such:
```julia
using Pkg
Pkg.add([
    "PlotlyJS",
    "CSV",
    "DataFrames",
    "DifferentialEquations",
    "Statistics",
    "HypothesisTests",
    "Distributions",
    "Random",
    "Base.Threads",
    "OrdinaryDiffEq",
    "StatsBase",
    "PosteriorStats"
])
```

## Workflow
### 1. Structural identifiability
The structural identifiability analysis is done with the MATLAB toolbox STRIKE-GOLDD. Install STRIKE-GOLDD from the software's [repository](https://github.com/afvillaverde/strike-goldd). Follow the instructions from the [user manual](https://github.com/afvillaverde/strike-goldd/blob/master/STRIKE-GOLDD/doc/STRIKE-GOLDD_manual.pdf), replacing the `models` file by the one from this repository.

## Citations
If you use any of the data or this code, please cite the associated publication. 

If you use STRIKE-GOLDD in your research, please cite the following paper:
Villaverde AF, Barreiro A, Papachristodoulou A. Structural identifiability of dynamic systems biology models. PLOS Computational Biology. 2016;12(10):e1005153. doi: 10.1371/journal.pcbi.1005153.
