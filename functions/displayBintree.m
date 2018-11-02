function displayBintree(bintree)
nBinNodes = size(bintree,1);
AM_bin_true = treeAMdirect(bintree);
labelsTrue = cell(1,nBinNodes);
for i = 1:nBinNodes
    labelsTrue{i} = ['v' int2str(i) '_l=' int2str(bintree(i,4)) '_f=' num2str(bintree(i,7)) '_t=' int2str(bintree(i,6)) '_o=' num2str(bintree(i,5))];
end
bg = biograph(AM_bin_true,labelsTrue);
view(bg)