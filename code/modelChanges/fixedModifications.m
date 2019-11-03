function model = fixedModifications(model,chemostat)
  % fixedModifications
  %   Additional modifications for improved simulations. Taken from GECKO:
  %   https://github.com/SysBioChalmers/GECKO/blob/2e5874714104722ed789da9caa333639b4b63069/Matlab_Module/limit_proteins/fixedModifications.m
  %
  %   model        (struct) metabolic model
  %   chemostat    (vector) acetate and pyruvate production rates
  %
  %   model        (struct) modified metabolic model
  %
  %   Usage: model = fixedModifications(model,chemostat)
  %

% Current values (aerobic Yeast 7.6)
NGAM = 0.7;

% Add NGAM reaction:
%            ATP  +  H2O  ->  ADP  +   H+   +  PO4
mets   = {'s_0434','s_0803','s_0394','s_0794','s_1322'};
coeffs = [-1,-1,1,1,1];
model  = addReaction(model,'NGAM', ...
    'metaboliteList', mets, ...
    'stoichCoeffList', coeffs, ...
    'reversible', false, ...
    'lowerBound', NGAM, ...
    'upperBound', NGAM);

if nargin == 2
    % Limit measured excretions:
    model.ub(strcmp(model.rxnNames,'acetate exchange'))  = chemostat(1);
    model.ub(strcmp(model.rxnNames,'pyruvate exchange')) = chemostat(2);

    % Limit unmeasured excretions:
    model.ub(strcmp(model.rxnNames,'(R,R)-2,3-butanediol exchange')) = 1e-5;
    model.ub(strcmp(model.rxnNames,'acetaldehyde exchange'))         = 1e-5;
end

end
