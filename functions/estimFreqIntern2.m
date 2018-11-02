function [freqDistrs,freqDistrLast,probMutEvents] = estimFreqIntern2(bintree,eps)
% bintree: [leftChild righChild parent label frequency timet fitness]
intern = find(bintree(:,1) + bintree(:,2) > 0);
nIntern = size(intern,1);
nNodes = size(bintree,1);
aux = [intern bintree(intern,6)];
aux = sortrows(aux,2);
internSorted = aux(:,1);
timesSorted = aux(:,2);
freqDistrs = zeros(nNodes,nNodes);
probMutEvents = zeros(1,nNodes);
probMutEvents(internSorted(1)) = 1;
freqDistrs(1,bintree(1,1)) = eps;
freqDistrs(1,bintree(1,2)) = 1-eps;

for i = 2:nIntern
    u = internSorted(i);
    p = internSorted(i-1);
    [freqDistrs(u,:),probMutEvents(u)] = estimFreqInternVert(bintree,u,p,freqDistrs(p,:),eps);
end

leafs = (find(bintree(:,1) + bintree(:,2) == 0))';
t = bintree(leafs(1),6) - timesSorted(end);
freqDistrLast = zeros(1,nNodes);
fm = max(bintree(leafs,7));
for j=leafs
    f = bintree(j,7);
    freqDistrLast(j) = freqDistrs(internSorted(end),j)*exp((f-fm)*t);
end
freqDistrLast = freqDistrLast/sum(freqDistrLast,2);
