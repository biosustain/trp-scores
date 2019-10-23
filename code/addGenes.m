
function geneTable = addGenes(geneTable,k_pos,new_set)

N = length(geneTable(:,1));
for i = 1:length(new_set.genes)
    pos = strcmp(geneTable(:,1),new_set.genes{i});
    if sum(pos) > 0
        geneTable{pos,k_pos} = new_set.k_genes(i);
    else
        N = N + 1;
        geneTable{N,1}     = new_set.genes{i};
        geneTable{N,k_pos} = new_set.k_genes(i);
    end
end

end
