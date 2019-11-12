function writePathwayCoverage(model,targets)
  % writePathwayCoverage
  %   Writes in ./results/pathwayCoverage.tsv a table with the coverage in each
  %   pathway of targets. Additionally, plots results for PPP and glycolysis.
  %
  %   model      (struct) metabolic model (in RAVEN format)
  %   targets    (cell) gene names that were detected as targets by compareDist.m
  %
  %   The output table has the following format:
  %   * 1st column: pathway name
  %   * 2nd column: number of target genes in pathway
  %   * 3rd column: number of total genes in pathway
  %   * 4th column: % of coverage in pathway
  %   * 5th column: significance of coverage (p-val)
  %
  %   Usage: writePathwayCoverage(model,targets)
  %

% Construct vector of pathways and vector of counts:
pathways = cell(1000,1);
counts   = zeros(1000,2);
x        = 1;
for i = 1:length(model.genes)
    % Find pathways associated to the gene:
    rxn_pos     = model.rxnGeneMat(:,i) ~= 0;
    pathway_set = model.subSystems(rxn_pos);
    pathway_set = cat(2, pathway_set{:});
    pathway_set = unique(pathway_set);

    for j = 1:length(pathway_set)
        if ~isempty(pathway_set{j})
            % Find pathway name:
            sce_pos = strfind(pathway_set{j},'sce');
            if isempty(sce_pos)
                pathway_name = pathway_set{j};
            else
                pathway_name = pathway_set{j}(sce_pos+10:end);
            end

            % Add new pathways the the array:
            if sum(strcmp(pathways,pathway_name)) == 0
                pathways{x,1} = pathway_name;
                x = x + 1;
            end
            pathway_pos = strcmp(pathways,pathway_name);

            % Count gene:
            if ismember(model.genes{i},targets)
                counts(pathway_pos,1) = counts(pathway_pos,1) + 1;
            end
            counts(pathway_pos,2) = counts(pathway_pos,2) + 1;
        end
    end
end

% Remove empty rows
pathways(x:end) = [];
counts(x:end,:) = [];

% Compute coverage %s & p-values:
counts(:,3) = counts(:,1)./counts(:,2)*100;
for i = 1:length(pathways)
    in_target   = counts(i,1);
    in_not      = counts(i,2) - in_target;
    out_target  = length(targets) - in_target;
    out_not     = length(model.genes) - length(targets) - in_not;
    Fmatrix     = [in_target,in_not;out_target,out_not];
    [~,p.val]   = fishertest(Fmatrix);
    counts(i,4) = p.val;
end

% Sort by p-values:
[counts,order] = sortrows(counts,4,'ascend');
pathways       = pathways(order);

% Write results:
fid = fopen('../../results/pathwayCoverage.tsv','wt');
fprintf(fid,'pathway\ttarget.count\tpathway.size\tcoverage.perc\tp.val\n');
for i = 1:length(pathways)
    fprintf(fid,'%s\t%d\t%d\t%.1f\t%.2e\n', pathways{i}, counts(i,:));
end
fclose(fid);

% Data to plot:
labels    = {'PPP','Glycolysis','Whole metabolism'};
pos_PPP   = strcmp(pathways,'Pentose phosphate pathway');
pos_glyco = strcmp(pathways,'Glycolysis');
data(1)   = counts(pos_PPP,3);
data(2)   = counts(pos_glyco,3);
data(3)   = length(targets)/length(model.genes)*100;

% Plot:
figure('position',[100,100,800,400])
hold on
b = bar(data);
b.FaceColor = 'flat';
b.CData(1,:) = [255,153,0]/255;
b.CData(2,:) = [0,32,76]/255;
b.CData(3,:) = [204,204,204]/255;
xlim([0.2,3.8])
ylim([0,100])
text_size = 15;
ylabel('Targets present in pathway [%]','FontSize',text_size);
set(gca,'XTick',1:3)
set(gca,'XTickLabel',labels)
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

end
