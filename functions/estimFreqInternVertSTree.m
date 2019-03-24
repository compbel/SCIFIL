function [freqDistr,probMutEvent,phi] = estimFreqInternVertSTree(u,h1,h2,t,fit,freqDistr,eps,usephi)
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
% is mutation 0
% h1 = stree{h1+1,2};
% h2 = stree{h2+1,2};
% n = length(freqDistr);
% freqDistr = zeros(1,n);
% currNodes = find(freqDistrPrev>0);
% fm = max(fit(currNodes));
fm = max(fit);
% for j=currNodes
% exps = exp(t*(fit - fm));
% freqDistr = freqDistrPrev.*exps;
totPop = 0;
for j=1:length(freqDistr)
    if freqDistr(j) > 0
%         f = fit(j);
        freqDistr(j) = freqDistr(j)*exp((fit(j)-fm)*t);
        totPop = totPop + freqDistr(j);
    end
end
% freqDistr = freqDistr/sum(freqDistr,2);
freqDistr = freqDistr/totPop;
probMutEvent = freqDistr(h2);
if usephi
    phi = freqDistr*fit';
else
    phi = 0;
end
freqDistr(h1) = eps*freqDistr(h2);
freqDistr(h2) = (1-eps)*freqDistr(h2);


