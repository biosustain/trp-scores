function [list,deleted] = deleteRepeated(list)
  % deleteRepeated
  %   Removes all repeated elements of a list.
  %
  %   list       (cell) a list
  %
  %   list       (cell) list with only non repeated elements
  %   deleted    (cell) list of the positions deleted
  %
  %   Usage: [list,deleted] = deleteRepeated(list)
  %

N       = length(list);
deleted = false(1,N);
for i = 1:N-1
    for j = i+1:N
        if strcmp(list{i},list{j})
            deleted(j) = true;
        end
    end
end
list(deleted) = [];

end
