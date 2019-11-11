
function writeGeneLabels(model,targets)
  % writeGeneLabels
  %   writes up a table according to the format:
  %   gene - pathway - is_target
  %   gene:         name of the gene
  %   pathway:      KEGG pathway the gene belongs to, based on the reactions the
  %                 gene is associated to. In case of more than one, the majority wins.
  %   is_target:    is it part of the target group (true/false)
  %
  %   model      (struct) metabolic model (in RAVEN format)
  %   targets    (cell) gene names that were detected as targets by compareDist.m
  %
  %   Usage: writeGeneLabels(model,targets)
  %

disp('Creating labels...')
fid = fopen('geneLabels.sif','wt');
fprintf(fid,'gene\tpathway\tis_target\n');

for i = 1:length(model.genes)
    gene = model.genes{i};

    % Determine pathway associated to gene:
    rxns     = model.rxnGeneMat(:,i) ~= 0;
    pathways = model.subSystems(rxns);
    is_glyco = findTarget(pathways,'sce00010');
    is_ppp   = findTarget(pathways,'sce00030');
    if is_ppp
        pathway = 'ppp';
    elseif is_glyco
        pathway = 'glycolysis';
    else
        pathway = 'other';
    end

    % Determine if gene is target or not:
    if ismember(model.genes{i},targets)
        is_target = 'true';
    else
        is_target = 'false';
    end

    % Write data:
    fprintf(fid,'%s\t%s\t%s\n', gene, pathway, is_target);
end

fclose(fid);

end

%%

function is_present = findTarget(pathways,target)

is_present = false;
for i = 1:length(pathways)
    for j = 1:length(pathways{i})
        if startsWith(pathways{i}{j},target)
            is_present = true;
        end
    end
end

end
