function plotPathways(model,targets)
  % plotPathways
  %   Plots a bar plot with the top 10 pathways with most targets
  %
  %   model      (struct) metabolic model (in RAVEN format)
  %   targets    (cell) gene names that were detected as targets by compareDist.m
  %
  %   Usage: plotPathways(model)
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

% Data to plot:
labels    = {'PPP','Glycolysis','Whole metabolism'};
pos_PPP   = strcmp(pathways,'Pentose phosphate pathway');
pos_glyco = strcmp(pathways,'Glycolysis');
data(1)   = counts(pos_PPP,1)/counts(pos_PPP,2)*100;
data(2)   = counts(pos_glyco,1)/counts(pos_glyco,2)*100;
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
