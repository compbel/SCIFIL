function [freqDistr,probMutEvent] = estimFreqInternVert(bintree,u,p,freqDistrPrev,eps)
% bintree: [leftChild righChild parent label frequency timet fitness]
nNodes = size(bintree,1);
freqDistr = zeros(1,nNodes);
currNodes = find(freqDistrPrev>0);
t = abs(bintree(u,6)-bintree(p,6));
%     freqDistrVPA = vpa(zeros(1,length(currNodes)),100);
fm = max(bintree(currNodes,7));
for j=currNodes
    f = bintree(j,7);
    freqDistr(j) = freqDistrPrev(j)*exp((f-fm)*t);
    if isinf(freqDistr(j))
        ['stop'];
    end
end
freqDistr = freqDistr/sum(freqDistr,2);
freqDistr(bintree(u,1)) = eps*freqDistr(u);
freqDistr(bintree(u,2)) = (1-eps)*freqDistr(u);
probMutEvent = freqDistr(u);
freqDistr(u) = 0;


