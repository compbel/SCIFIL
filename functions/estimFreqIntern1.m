function [freqDistrs,freqDistrLast] = estimFreqIntern1(bintree,eps)
% bintree: [leftChild righChild parent label frequency timet fitness]
intern = find(bintree(:,1) + bintree(:,2) > 0);
nIntern = size(intern,1);
nNodes = size(bintree,1);
aux = [intern bintree(intern,6)];
aux = sortrows(aux,2);
internSorted = aux(:,1);
timesSorted = aux(:,2);
freqDistrs = zeros(nNodes,nNodes);
freqDistrs(1,bintree(1,1)) = eps;
freqDistrs(1,bintree(1,2)) = 1-eps;
nodeBirth = zeros(1,nNodes);
for i = 2:nNodes
    nodeBirth(i) = bintree(bintree(i,3),6);
end

for i = 2:nIntern
    u = internSorted(i);
    p = internSorted(i-1);
    currNodes = intersect(find(nodeBirth <= bintree(u,6)),find(freqDistrs(p,:)>0));
    t = timesSorted(i) - timesSorted(i-1);
%     freqDistrVPA = vpa(zeros(1,length(currNodes)),100);
    fm = max(bintree(currNodes,7));
    for j=currNodes
        f = bintree(j,7);
        freqDistrs(u,j) = freqDistrs(p,j)*exp((f-fm)*t);
        if isinf(freqDistrs(u,j))
            ['stop'];
        end
    end
    freqDistrs(u,:) = freqDistrs(u,:)/sum(freqDistrs(u,:),2);
    freqDistrs(u,bintree(u,1)) = eps*freqDistrs(u,u);
    freqDistrs(u,bintree(u,2)) = (1-eps)*freqDistrs(u,u);
%     freqDistrs(u,currNodes) = calcFreqDistInterv(freqDistrs(p,currNodes)',bintree(currNodes,7),t);
%     freqDistrs(u,bintree(u,1)) = eps;
%     freqDistrs(u,bintree(u,2)) = max(0,freqDistrs(u,u)-eps);
%     freqDistrs(u,bintree(u,2)) = freqDistrs(u,u)-eps;
    freqDistrs(u,u) = 0;
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
