function [likelihood,fit] = probTreeParamOrderSTreePartial1(stree,internOrder,theta,obsFreqLeafs,eps)
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
% is mutation 0

nIntern = length(internOrder);
wildType = stree{1,2};
nLeafs = size(stree,1) - nIntern -1;

fit = zeros(nLeafs,1);
fit(wildType) = 1;
for k = 1:length(internOrder)
    u = internOrder(k) +1; % +1 ???
    p = stree{u,3};
    i = stree{u,2};
    j = stree{p,2};
    o1 = obsFreqLeafs(i);
    o2 = obsFreqLeafs(j);
    fit(i) = fit(j) - theta*(log(eps/(1-eps)) + log(o2/o1))/(nIntern-k+1);
end

fit = fit';

T = 10^7;
stree = estimTimeSTree(stree,fit,obsFreqLeafs,internOrder,T,eps);
[freqDistInt,freqDistLeaf,probMutEvents] = estimFreqInternSTree(stree,fit,eps);
% Height = max([stree{:,6}]);


likelihood = 0;
for i = 2:length(internOrder)
    u = internOrder(i);
    v = internOrder(i-1);
%     freq = freqDistInt(u,x) + freqDistInt(u,y);
    freq = probMutEvents(u);
    t = stree{u+1,6} - stree{v+1,6};
%     pois = poisspdf(round(t),1/theta);
%     likelihood = likelihood + log(freq) + log(pois);
    likelihood = likelihood + log(freq);
end
% t = Height - stree{internOrder(end)+1,6};
% pois = poisspdf(round(t),1/theta);
% likelihood = likelihood + log(pois);
