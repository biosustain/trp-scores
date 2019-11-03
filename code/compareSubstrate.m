function FC = compareSubstrate(model,product,substrate)
  % compareSubstrate
  %   Performs a full FSEOF analysis for a specific product/substrate pair.
  %
  %   model        (struct) metabolic model (in irreversible format)
  %   product      (str) name of exchange reaction to optimize for
  %   substrate    (str) either "glucose" or "ethanol"
  %
  %   FC           (struct) all results, including:
  %                * FC.flux_WT: Flux distribution for the WT strain (100% of carbon towards growth).
  %                * FC.rxns: A list with all reactions with fluxes that change
  %                           consistently as tryptophan production increases.
  %                           Column format: ID - name - gene association - chemical equation
  %                * FC.k_rxns: Scores for each of the reactions in "rxns"
  %                * FC.v_matrix: Fluxes for each reaction in "rxns" and each of the 11 simulated conditions.
  %                * FC.k_matrix: Fold changes for each reaction in "rxns" and each of the 11 simulated conditions.
  %                * FC.genes: Names of the selected targets.
  %                * FC.k_genes: Scores for each of the targets in "genes".
  %
  %   Usage: FC = compareSubstrate(model,product,substrate)
  %

% Simulate WT (100% growth) and forced (X% growth and the rest towards product):
FC.flux_WT = simulateGrowth(model,product,substrate,1);
alpha      = 0.3:0.05:0.8;
v_matrix   = zeros(length(model.rxns),length(alpha));
k_matrix   = zeros(length(model.rxns),length(alpha));
for i = 1:length(alpha)
    flux_MAX      = simulateGrowth(model,product,substrate,alpha(i));
    v_matrix(:,i) = flux_MAX;
    k_matrix(:,i) = flux_MAX./FC.flux_WT;
end

% Generate rxn equations:
rxnEqs = printRxnFormula(model,model.rxns,true,true,true);

% Take out rxns with no grRule:
withGR   = ~cellfun(@isempty,model.grRules);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
gene_rxn = model.rxnGeneMat(withGR,:);
FC.rxns  = [model.rxns(withGR) model.rxnNames(withGR) model.grRules(withGR) rxnEqs(withGR)];

% Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
gene_rxn = gene_rxn(non_nan,:);
FC.rxns  = FC.rxns(non_nan,:);

% Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

% Replace any Inf value with 1000 (maximum value is ~700):
k_matrix(isinf(k_matrix)) = 1000;

% Filter out values that are inconsistent at different alphas:
always_down  = sum(k_matrix <= 1,2) == length(alpha);
always_up    = sum(k_matrix >= 1,2) == length(alpha);
incons_rxns  = always_down + always_up == 0;
incons_genes = sum(gene_rxn(incons_rxns,:),1) > 0;
incons_rxns  = sum(gene_rxn(:,incons_genes),2) > 0;
v_matrix     = v_matrix(~incons_rxns,:);
k_matrix     = k_matrix(~incons_rxns,:);
gene_rxn     = gene_rxn(~incons_rxns,:);
FC.rxns      = FC.rxns(~incons_rxns,:);

% Order from highest to lowest k:
FC.k_rxns   = mean(k_matrix,2);
[~,order]   = sort(FC.k_rxns,'descend');
FC.k_rxns   = FC.k_rxns(order,:);
FC.v_matrix = v_matrix(order,:);
FC.k_matrix = k_matrix(order,:);
gene_rxn    = gene_rxn(order,:);
FC.rxns     = FC.rxns(order,:);

% Create list of remaining genes and filter out any inconsistent score:
FC.genes   = model.genes(sum(gene_rxn,1) > 0);
FC.k_genes = zeros(size(FC.genes));
gene_rxn   = gene_rxn(:,sum(gene_rxn,1) > 0);
cons_genes = false(size(FC.genes));
for i = 1:length(FC.genes)
    k_set         = FC.k_rxns(gene_rxn(:,i) > 0);
    always_down   = sum(k_set <= 1) == length(k_set);
    always_up     = sum(k_set >= 1) == length(k_set);
    cons_genes(i) = always_down + always_up == 1;
    FC.k_genes(i) = mean(k_set);
end
FC.genes   = FC.genes(cons_genes);
FC.k_genes = FC.k_genes(cons_genes);

% Filter any value between mean(alpha) and 1:
unchanged  = (FC.k_genes >= mean(alpha) - 1e-3) + (FC.k_genes <= 1 + 1e-3) == 2;
FC.genes   = FC.genes(~unchanged);
FC.k_genes = FC.k_genes(~unchanged);

% Order from highest to lowest k:
[~,order]  = sort(FC.k_genes,'descend');
FC.genes   = FC.genes(order,:);
FC.k_genes = FC.k_genes(order,:);

end
