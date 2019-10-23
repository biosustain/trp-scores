
%Load COBRA and models:
initCobraToolbox
model = load('../model/yeast_7.6_cobra.mat');
model = model.model;

%Make basic changes to model:
cd modelChanges
model = changeModel(model);
cd ..

%Block tryptophan reutilization (BNA2):
model.ub(strcmp(model.rxns,'r_0694')) = 0; %L-trp + O2 -> N-formyl-L-kynurenine

%Block alternative phosphofruktokinase:
model.ub(strcmp(model.rxns,'r_0887')) = 0; %ATP + sedohept-7P -> ADP + H+ + sedohept-1,7biP

results.model    = model;
results.glucose  = compareSubstrate(model,'L-tryptophan exchange','glucose');
results.ethanol  = compareSubstrate(model,'L-tryptophan exchange','ethanol');

%Create gene table:
results.geneTable      = cell(length(results.glucose.genes),3);
results.geneTable(:,1) = results.glucose.genes;
results.geneTable(:,2) = num2cell(results.glucose.k_genes);
results.geneTable      = addGenes(results.geneTable,3,results.ethanol);
