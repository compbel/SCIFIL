function AM_tree = treeAMdirect(tree)
n = size(tree,1);
AM_tree = zeros(n,n);
for i = 1:n
    child1 = tree(i,1);
    child2 = tree(i,2);
    if child1 > 0
        AM_tree(i,child1) = 1;
    end
    if child2 > 0
        AM_tree(i,child2) = 1;
    end
end