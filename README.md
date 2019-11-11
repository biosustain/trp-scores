
# trp-scores

Code used for scoring genes as upregulation / downregulation candidates for increased tryptophan production.

### Requirements

* A functional Matlab installation (MATLAB 7.3 or higher).
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).

### Usage

* Run `./code/compareDist.m` for performing the complete target search.
* Run `./code/networkVisualization/networkVisualization.m` (using `results` as input) to generate the required `.sif` files for constructing the associated cytoscape files.

### Documentation

In the Matlab console, type `help function_name` for a description of each Matlab function.
