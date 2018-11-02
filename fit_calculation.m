load('data_for_validation_with_noise.mat');

err = zeros(1,length(data));
correlation = zeros(1, length(data));
l_infer = zeros(1,length(data));
l_real = zeros(1,length(data));
global gl_oter;
global l_min;
global time_limit;
gl_oter = 1;
time_limit = 10 * 60;
for iter = 1:length(data)
    l_min = 10e8;
    gk_oter = 1;
    bintree = data(iter).bintree;
    bintreeTrue = data(iter).bintreeTrue;
    cellDistr = data(iter).cellDistr;
    M = data(iter).M;
    leaf_count = 0;
    for i = 1:size(bintreeTrue, 1)
        if bintreeTrue(i,1)+bintreeTrue(i,2) == 0
            leaf_count = leaf_count +1;
        end
    end
    obsFreqLeafs = cellDistr/sum(cellDistr);
    l_index = 1;
    real_fit = zeros(1,leaf_count);
    for i = 1:size(bintreeTrue, 1)
        if bintreeTrue(i,1)+bintreeTrue(i,2) == 0
            real_fit(l_index)= bintreeTrue(i,7);
            l_index = l_index + 1;
        end
    end
    
%     [likelihood,fit] = bintreeFit(bintree,theta,f,obsFreqLeafs,eps,Tmax,M,M,pFalseNeg,pFalsePos);
%      [likelihood1,fit1, dist_table] = bintreeFitNew(bintree,theta,f,obsFreqLeafs,eps,Tmax);
    [likelihood,fit, order,stats] = bintreeFitBruteForce(bintree,theta,f,obsFreqLeafs,eps,Tmax, time_limit);
    nIntern = size(M,2);
    nUnCells = size(M,1);
    stree = bintreeToStree(bintree, nIntern,nUnCells, cellDistr);
    disp('start stree')
    [likelihoods,fits, orders,statss] = fitBruteForceStree(stree,theta,f,obsFreqLeafs,eps,Tmax);
    di = likelihoods-likelihood
    l_infer(iter) = likelihood;
    l_real(iter) = -probTreeParam1(bintree,theta,real_fit,obsFreqLeafs,eps,Tmax, 0);
    err(iter) = sqrt(immse(fit, real_fit));
    correlation(iter) = corr(fit.',real_fit.', 'type','Spearman');
    disp(iter);
end

