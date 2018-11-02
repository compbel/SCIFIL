function [likelihood, fit, order,stats] = bintreeFitBruteForce(bintree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax, timeLimit)

global gl_oter;
global l_min;
global bestOrder;
global filterFit;
global filterLB;
global time_limit;
l_min = 10e8;
gl_oter = 1;
bestOrder = [];
filterFit = 0;
filterLB = 0;
time_limit = timeLimit;

    leafs = (find(bintree(:,1)+bintree(:,2) == 0))';
    leaf_count = 0;
    for i = 1:size(bintree, 1)
        if bintree(i,1)+bintree(i,2) == 0
            leaf_count = leaf_count +1;
            bintree(i, 5) = obsFreqLeafs(leaf_count);
        end
    end
    zeromut = 1;
    while bintree(zeromut,2) ~= 0
        zeromut = bintree(zeromut, 2);
    end
    bintree(zeromut, 7) = 1;
    nodes_count = size(bintree, 1) - leaf_count;
    [leftHapl,rightHapl] = getProducedHapl(bintree);
    freqDistr = zeros(1, nodes_count+leaf_count);
    
    freqDistr(bintree(1,2)) = 1-eps;
    freqDistr(bintree(1,1)) = eps;
    tic;
    [likelihood, fit] = recEdgeLengthEstimation(1, bintree, [], [1], theta, obsFreqLeafs,eps,Tmax, leftHapl, rightHapl, nodes_count, fmax, fmin, freqDistr, 1);
    toc
    order = bestOrder;
    stats = [gl_oter filterFit filterLB];
