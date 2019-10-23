
# trp-scores

Code used for scoring genes as upregulation / downregulation candidates for increased tryptophan production.

### Requirements

* A functional Matlab installation (MATLAB 7.3 or higher).
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).

### Usage

Run `./code/compareDist.m` for performing the complete analysis. The `results` variable will contain all relevant information from the analysis, including:

* `model`: The modified model.
* `glucose`: Results growing on glucose.
* `ethanol`: Results growing on ethanol.
* `geneTable`: Merged gene table with both glucose and ethanol results.

Both `glucose` and `ethanol` wil in turn contain:

* `flux_WT`: Flux distribution for the WT strain (100% of carbon towards growth).
* `rxns`: A list with all reactions with fluxes that change consistently as tryptophan production increases. Column format: ID - name - gene association - chemical equation
* `k_rxns`: Scores for each of the reactions in `rxns`
* `v_matrix`: Fluxes for each reaction in `rxns` and each of the 11 simulated conditions.
* `k_matrix`: Fold changes for each reaction in `rxns` and each of the 11 simulated conditions.
* `genes`: Names of the selected targets.
* `k_genes`: Scores for each of the targets in `genes`.
