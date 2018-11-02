function [likelihood, fit, order,stats] = fitBruteForceStree(stree,theta,fmax,obsFreqLeafs,eps,Tmax)

global gl_oter;
global l_min;
global bestOrder;
global filterFit;
global filterLB;
l_min = 10e8;
gl_oter = 1;
bestOrder = [];
filterFit = 0;
filterLB = 0;


leaf_count = 0;
for i = 1:size(stree, 1)
    if isempty(stree{i,2})
        leaf_count = leaf_count +1;
        stree{i, 5} = obsFreqLeafs(leaf_count);
    end
end
nodes_count = size(stree, 1) - leaf_count - 1;
stree{1,7} = 1;
stree{1,6} = 0;
tic
[likelihood, fit] = recEdgeLengthEstimationStree(1, stree, [], stree{1,1}, [], theta, obsFreqLeafs,eps,Tmax, nodes_count, fmax, [], 1);
toc
order = bestOrder;
stats = [gl_oter filterFit filterLB];
