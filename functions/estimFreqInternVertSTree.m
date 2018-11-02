function [freqDistr,probMutEvent] = estimFreqInternVertSTree(u,h1,h2,t,fit,freqDistrPrev,eps)
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
% is mutation 0
freqDistr = zeros(1,length(freqDistrPrev));
currNodes = find(freqDistrPrev>0);
fm = max(fit(currNodes));
for j=currNodes
    f = fit(j);
    freqDistr(j) = freqDistrPrev(j)*exp((f-fm)*t);
    if isinf(freqDistr(j))
        ['stop'];
    end
end
freqDistr = freqDistr/sum(freqDistr,2);
probMutEvent = freqDistr(h2);
freqDistr(h1) = eps*freqDistr(h2);
freqDistr(h2) = (1-eps)*freqDistr(h2);



