function networkVisualization(results)
  % networkVisualization
  % Generates visualization files.
  %
  %   model      (struct) metabolic model (in COBRA format)
  %   results    (cell) results from compareDist.m
  %
  %   Usage: networkVisualization(results)
  %

% Pre-process input:
targets = results.geneTable(:,1);
model   = ravenCobraWrapper(results.model);
format short

% Load newest version of the yeast model:
outfilename = websave('yeastGEM.mat','https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/mat/yeastGEM.mat');
model_new   = load(outfilename);
model_new   = ravenCobraWrapper(model_new.model);
delete(outfilename)

% Assign subsystems to old model:
for i = 1:length(model.rxns)
    pos_new = strcmp(model_new.rxns,model.rxns{i});
    if sum(pos_new) > 0
        model.subSystems(i) = model_new.subSystems(pos_new);
    end
end

% Write pathway coverage file & plot results for PPP + glycolysis:
writePathwayCoverage(model,targets)

% Write gene labels (indicating pathway and target):
writeGeneLabels(model,targets)

% Remove currency metabolites:
currency_mets = {'H2O','H+','carbon dioxide','oxygen','phosphate', ...
                 'diphosphate','ammonium','ATP','ADP','AMP','NAD', ...
                 'NAD(+)','NADH','NADP(+)','NADPH'};
model = removeMets(model,currency_mets,true);

% Write gene interactions based on shared metabolites between genes:
writeGeneInteractions(model)

end
