
# trp-scores

Code for scoring genes as upregulation / downregulation candidates for increased tryptophan production. Used in [this publication](https://doi.org/10.1038/s41467-020-17910-1).

### Requirements

* A functional Matlab installation (MATLAB 7.3 or higher).
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).

### Usage

* Run `./code/compareDist.m` for performing the complete target search.
* Run `./code/networkVisualization/networkVisualization.m` (using `results` as input) to generate the required `.sif` files for constructing the associated cytoscape files.

### Documentation

In the Matlab console, type `help function_name` for a description of each Matlab function.

### Citation

Jie Zhang, Søren D. Petersen, Tijana Radivojevic, Andrés Ramirez, Andrés Pérez-Manríquez, Eduardo Abeliuk, Benjamín J. Sánchez, Zak Costello, Yu Chen, Michael J. Fero, Hector Garcia Martin, Jens Nielsen, Jay D. Keasling & Michael K. Jensen (2020). _Combining mechanistic and machine learning models for predictive engineering and optimization of tryptophan metabolism._ Nat Commun 11, 4880; doi: https://doi.org/10.1038/s41467-020-17910-1

