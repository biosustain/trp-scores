
function writeGeneInteractions(model)
  % writeGeneInteractions
  %   Creates gene-gene interaction file in /results/geneInteractions.sif, which
  %   will have a connection size equal to the number of metabolites shared by
  %   the associated reactions of the 2 genes.
  %
  %   model    (struct) metabolic model (in RAVEN format)
  %
  %   Usage: writeGeneInteractions(model)
  %

fid = fopen('../../results/geneInteractions.sif','wt');
for i = 1:length(model.genes)-1
    for j = i+1:length(model.genes)
        % Genes:
        gene_i = model.genes{i};
        gene_j = model.genes{j};

        % Reactions each gene is associated to:
        rxns_i = model.rxnGeneMat(:,i) ~= 0;
        rxns_j = model.rxnGeneMat(:,j) ~= 0;

        % Metabolites those reactions are associated to:
        mets_i = sum(model.S(:,rxns_i) ~= 0, 2);
        mets_j = sum(model.S(:,rxns_j) ~= 0, 2);

        % Shared metabolites:
        shared_mets = full(sum((mets_i > 0).*(mets_j > 0)));
        if shared_mets > 0
            fprintf(fid,'%s\t%.1f\t%s\n', gene_i, shared_mets, gene_j);
        end
    end
    if mod(i,100) == 0
      disp(['Creating interaction file: ready with gene #' num2str(i)])
    end
end
fclose(fid);

end
