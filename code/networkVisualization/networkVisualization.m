function networkVisualization(results)
  % networkVisualization
  % generates .sif files for visualization of the network in cytoscape
  %
  %   model      (struct) metabolic model (in COBRA format)
  %   results    (cell) results from compareDist.m
  %
  %   Usage: networkVisualization(model,results)
  %

% Pre-process input:
targets   = results.geneTable(:,1);
model_old = ravenCobraWrapper(results.model);
format short

% Load newest version of the yeast model:
outfilename = websave('yeastGEM.mat','https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/mat/yeastGEM.mat');
model_new   = load(outfilename);
model_new   = ravenCobraWrapper(model_new.model);
delete(outfilename)

% Write gene labels using the latest yeast-GEM model (as the old one does not
% have pathway information):
writeGeneLabels(model_new,targets)

% Write gene interactions in the model with which the analysis was made:
% TODO: test removing currency metabolites
writeGeneInteractions(model_old)

end
