
# trp-scores

Code for scoring genes as upregulation / downregulation candidates for increased tryptophan production. Used in [this pre-print](https://doi.org/10.1101/858464).

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

Jie Zhang, Søren D. Petersen, Tijana Radivojevic, Andrés Ramirez, Andrés Pérez, Eduardo Abeliuk, Benjamín J. Sánchez, Zachary Costello, Yu Chen, Mike Fero, Hector Garcia Martin, Jens Nielsen, Jay D. Keasling & Michael K. Jensen (2019). _Predictive engineering and optimization of tryptophan metabolism in yeast through a combination of mechanistic and machine learning models._ bioRxiv 858464; doi: https://doi.org/10.1101/858464
