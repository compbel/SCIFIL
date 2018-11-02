function [likelihood,fit] = probTreeParamOrder(bintree,internOrder,leftHapl,rightHapl,theta,obsFreqLeafs,eps)
% bintree: [leftChild righChild parent label frequency timet fitness]

nVert = size(bintree,1);
nLeafs = nVert - length(internOrder);
nIntern = length(internOrder);
wildType = rightHapl(1);
leafs = find(bintree(:,1) + bintree(:,2) == 0);

A = zeros(nLeafs,nLeafs);
r = zeros(nLeafs,1);
for k = 1:length(internOrder)
    u = internOrder(k);
    i = leftHapl(u);
    j = rightHapl(u);
    i1 = find(leafs == i);
    j1 = find(leafs == j);
    o1 = obsFreqLeafs(i1);
    o2 = obsFreqLeafs(j1);
    A(k,i1) = -1;
    r(k) = theta*(log(eps/(1-eps)) + log(o2/o1))/(nIntern-k+1);
    A(k,j1) = 1;
end
A(nLeafs,leafs == wildType) = 1;
r(nLeafs) = 1;

fit = linsolve(A,r);
fit = fit';

bintree = assignFit(bintree,obsFreqLeafs,fit);
T = 10^7;
bintreeInfer = estimTime1(bintree,T,eps); 
[freqDistInt,freqDistLeaf,probMutEvents] = estimFreqIntern2(bintreeInfer,eps);
Height = max(bintreeInfer(:,6));


likelihood = 0;
for i = 2:length(internOrder)
    u = internOrder(i);
    v = internOrder(i-1);
    x = bintreeInfer(u,1);
    y = bintreeInfer(u,2);
%     freq = freqDistInt(u,x) + freqDistInt(u,y);
    freq = probMutEvents(u);
    t = bintreeInfer(u,6) - bintreeInfer(v,6);
    pois = poisspdf(round(t),1/theta);
    likelihood = likelihood + log(freq) + log(pois);
end
t = Height - bintreeInfer(internOrder(end),6);
pois = poisspdf(round(t),1/theta);
likelihood = likelihood + log(pois);
