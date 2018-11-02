function [leftHapl,rightHapl] = getProducedHapl(bintree)
% bintree: [leftChild righChild parent label frequency timet fitness]
nVert = size(bintree,1);
leftHapl = zeros(1,nVert);
rightHapl = zeros(1,nVert);
leafs = find(bintree(:,1) + bintree(:,2) == 0);


AM = treeAMdirect(bintree);
G = digraph(AM);
sort = flip(toposort(G));

for v = sort
    if bintree(v,1) + bintree(v,2) == 0
        leftHapl(v) = v;
        rightHapl(v) = v;
    else 
        leftHapl(v) = rightHapl(bintree(v,1));
        rightHapl(v) = rightHapl(bintree(v,2));
    end
end
