function flux = simulateGrowth(model,rxn,substrate,alpha)
  % simulateGrowth
  %   Simulates the metabolic model for a suboptimal growth and all the rest of carbon
  %   towards the product of interest.
  %
  %   model        (struct) metabolic model (in irreversible format)
  %   rxn          (str) name of exchange reaction to optimize for
  %   substrate    (str) either "glucose" or "ethanol"
  %   alpha        (float) fraction between 0 and 1 representing the suboptimal growth
  %
  %   flux         (vector)
  %
  %   Usage: flux = simulateGrowth(model,rxn,substrate,alpha)
  %

% Block/allow consumptions:
if strcmp(substrate,'glucose')
    model.ub(strcmp(model.rxnNames,'D-glucose exchange'))              = 0;
    model.ub(strcmp(model.rxnNames,'D-glucose exchange (reversible)')) = 1;
    model.ub(strcmp(model.rxnNames,'ethanol exchange'))                = 1000;
    model.ub(strcmp(model.rxnNames,'ethanol exchange (reversible)'))   = 0;
elseif strcmp(substrate,'ethanol')
    model.ub(strcmp(model.rxnNames,'D-glucose exchange'))              = 1000;
    model.ub(strcmp(model.rxnNames,'D-glucose exchange (reversible)')) = 0;
    model.ub(strcmp(model.rxnNames,'ethanol exchange'))                = 0;
    model.ub(strcmp(model.rxnNames,'ethanol exchange (reversible)'))   = 1;
end

% Positions of biomass & target rxn:
posX = strcmp(model.rxnNames,'growth');
posP = strcmp(model.rxnNames,rxn);

% Max growth:
sol = optModel(model,posX,+1);

% Fix growth suboptimal and then max product:
model.lb(posX) = sol.x(posX)*0.999*alpha;
sol            = optModel(model,posP,+1);

% Fix also product and minimize fluxes:
model.lb(posP) = sol.x(posP)*0.999;
sol            = optModel(model,1:length(model.rxns),-1);
flux           = sol.x;

end

%%

function sol = optModel(model,pos,c,base_sol)
  % optModel
  %   Wrapper for optimizeCbModel with added functionality for irreversible models.
  %
  %   model       (struct) metabolic model (in irreversible format)
  %   pos         (int) the position of the reaction to optimize
  %   c           (int) either -1 (for minimizing) or +1 (for maximizing)
  %   base_sol    (struct, opt) a previous COBRA solution (for fixing values)
  %
  %   sol         (struct) COBRA solution from running the optimization
  %
  %   Usage: sol = optModel(model,pos,c,base_sol)
  %

% Make sure that the opposite reaction is blocked:
if length(pos) == 1
    rxn_code = model.rxns{pos};
    if strcmp(rxn_code(end-3:end),'_REV')
        rev_pos = strcmp(model.rxns,rxn_code(1:end-4));
    else
        rev_pos = strcmp(model.rxns,[rxn_code '_REV']);
    end
    model.lb(rev_pos) = 0;
    model.ub(rev_pos) = 0;
end

% Optimize:
model.c      = zeros(size(model.rxns));
model.c(pos) = c;
sol          = optimizeCbModel(model);

% If optimization didn't work, fix opposite reaction at a basal solution:
if isempty(sol.x) && length(pos) == 1
    model.lb(rev_pos) = base_sol.x(rev_pos);
    model.ub(rev_pos) = base_sol.x(rev_pos);
    sol               = optimizeCbModel(model);
end

end
