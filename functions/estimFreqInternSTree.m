function [freqDistrs,freqDistrLast,probMutEvents] = estimFreqInternSTree(stree,fit,eps)
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
% is mutation 0
intern = find(~cellfun(@isempty,stree(:,2)));
intern = intern(2:end);
leafs = find(cellfun(@isempty,stree(:,2)));
nIntern = size(intern,1);
nNodes = size(stree,1)-1;
nLeafs = size(stree,1) - nIntern - 1;

aux = [intern cell2mat(stree(intern,6))];
aux = sortrows(aux,2);
internSorted = aux(:,1);
timesSorted = aux(:,2);
freqDistrs = zeros(nIntern+1,nLeafs);
probMutEvents = zeros(1,nIntern);
probMutEvents(internSorted(1)-1) = 1;

u = internSorted(1)-1;
p = 0;

freqDistrs(u,stree{u+1,2}) = eps;
freqDistrs(u,stree{p+1,2}) = 1-eps;

for i = 2:nIntern
    u = internSorted(i)-1;
    p = internSorted(i-1)-1;
    h1 = stree{u+1,2};
    h2 = stree{stree{u+1,3},2};
    t = abs(stree{u+1,6}-stree{p+1,6});
    [freqDistrs(u,:),probMutEvents(u)] = estimFreqInternVertSTree(u,h1,h2,t,fit,freqDistrs(p,:),eps);
end

t = stree{leafs(1),6} - timesSorted(end);
freqDistrLast = zeros(1,nLeafs);
fm = max(fit);
for j=1:nLeafs
    f = fit(j);
    freqDistrLast(j) = freqDistrs(internSorted(end)-1,j)*exp((f-fm)*t);
end
freqDistrLast = freqDistrLast/sum(freqDistrLast,2);
