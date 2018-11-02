function bintree = assignFit(bintree,obsFreqLeafs,fit)
% bintree: [leftChild righChild parent label frequency timet fitness]
leafs = (find(bintree(:,1) + bintree(:,2) == 0))';
bintree(leafs,7) = fit;
bintree(leafs,5) = obsFreqLeafs;

AM = treeAMdirect(bintree);
G = digraph(AM);
sort = flip(toposort(G));

for v = sort
    if bintree(v,1) + bintree(v,2) ~= 0
        bintree(v,7) = bintree(bintree(v,2),7);
        bintree(v,5) = bintree(bintree(v,2),5);
    end
end