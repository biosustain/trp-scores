% compareDist
%   Main script of the analysis, for obtaining a list of genes that could
%   affect the production of tryptophan. After running it, the `results` variable
%   will contain all relevant information from the analysis, including:
%
%   * results.model: (struct) The modified model
%   * results.glucose: (struct) Results growing on glucose
%   * results.ethanol: (struct) Results growing on ethanol
%   * results.geneTable: (cell) Merged gene table with both glucose and ethanol results
%
%   A description of the "glucose" and "ethanol" structures is available at
%   compareSubstrate.m
%

% Load COBRA and models:
initCobraToolbox
model = load('../model/yeast_7.6_cobra.mat');
model = model.model;

% Make basic changes to model:
cd modelChanges
model = changeModel(model);
cd ..

% Block tryptophan reutilization (BNA2):
model.ub(strcmp(model.rxns,'r_0694')) = 0; %L-trp + O2 -> N-formyl-L-kynurenine

% Block alternative phosphofruktokinase:
model.ub(strcmp(model.rxns,'r_0887')) = 0; %ATP + sedohept-7P -> ADP + H+ + sedohept-1,7biP

results.model    = model;
results.glucose  = compareSubstrate(model,'L-tryptophan exchange','glucose');
results.ethanol  = compareSubstrate(model,'L-tryptophan exchange','ethanol');

% Create gene table:
results.geneTable      = cell(length(results.glucose.genes),3);
results.geneTable(:,1) = results.glucose.genes;
results.geneTable(:,2) = num2cell(results.glucose.k_genes);
results.geneTable      = addGenes(results.geneTable,3,results.ethanol);
