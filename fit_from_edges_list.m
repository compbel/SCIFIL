addpath('functions');
tree = convertEdgesListToTree(input);
converter = TreeConverter;
bintree = converter.convertToBinTree(tree);
fmax = 1.5;
fmin = 1;
theta = 0.01;
eps = 1e-5;
Tmax = 150000;
leafs = (find(bintree(:,1)+bintree(:,2) == 0))';
cellDist = bintree(leafs,5).';
obsFreqLeafs = cellDist/sum(cellDist);
timeLimit = 600;

% stree = bintreeToStree(bintree, size(bintree,1)-length(leafs),length(leafs), cellDist);
% intern = find(bintree(:,1) + bintree(:,2) > 0);
% orderTrue = intern;
% [leftHapl,rightHapl] = getProducedHapl(bintree);
% [likelihood2,fitInfer1] = probTreeParamOrder(bintree,orderTrue,leftHapl,rightHapl,theta,obsFreqLeafs,eps);
% [likelihoods,fits, orders,statss] = fitBruteForceStree(stree,theta,f,obsFreqLeafs,eps,Tmax);
[likelihood,fit,order,stats] = bintreeFitBruteForce(bintree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax,timeLimit);
% fill all fits according to initial order of mutations
true_fit_order = zeros(1,length(leafs)-1);
for i=1:length(leafs)
    if fit(i) == 1
        continue
    end
    % while current is a right child go up
    cur = leafs(i);
    par = bintree(cur, 3);
    while(bintree(par,2) == cur)
        cur = par;
        par = bintree(cur, 3);
    end
    % for parent node this node is left child so parent node mutation is
    % what we are looking for
    true_fit_order(bintree(par, 4)) = fit(i);
end
disp("")
disp("Found mutation order:");
disp(bintree(order,4).');
disp("Found fitness vector:");
disp(true_fit_order);

fileID = fopen(output,'w');
fprintf(fileID,'%u ', bintree(order,4).');
fprintf(fileID,'\n');
fprintf(fileID,'%6.4f ', true_fit_order);
fprintf(fileID,'\n');
fclose(fileID);
