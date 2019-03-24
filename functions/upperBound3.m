function [ubCut,ub,addOrder,fitLB] = upperBound3(stree,curFit,currLikelihood,chain1,chain2,lastCh1,lastCh2,lastOrder,freqCoeff,theta,freqDistr,bestLikelihood,eps)
% stree: [nextChild haplotype parent label frequency timet fitness oldChild]
ubCut = false;
ub = currLikelihood;
fitLB = curFit;
n = length(curFit) - 1;
freqDistrCurr = freqDistr;
addOrder = zeros(1,length(chain1)-lastCh1 + length(chain2) - lastCh2);
iChain1 = lastCh1+1;
iChain2 = lastCh2+1;
j = lastOrder;
count = 1;
while (iChain1 <= length(chain1)) && (iChain2 <= length(chain2))
    j = j + 1;
%     t = (j-1)/theta;
    p1 = stree{chain1(iChain1),3};
    f1 =  fitLB(p1) - freqCoeff(chain1(iChain1))*theta/(n-j+1);
    p2 = stree{chain2(iChain2),3};
    f2 =  fitLB(p2) - freqCoeff(chain2(iChain2))*theta/(n-j+1);
    if f1 < f2
        fitLB(chain1(iChain1)) = f1;
        [freqDistrCurr,probMutEvent,phi] = estimFreqInternVertSTree(0,stree{chain1(iChain1),2},stree{p1,2},1/theta, fitLB, freqDistrCurr,eps,false);
        p = p1;
        addOrder(count) = chain1(iChain1);
        count = count+1;
        iChain1 = iChain1 + 1;
    else
        fitLB(chain2(iChain2)) = f2;
        [freqDistrCurr,probMutEvent,phi] = estimFreqInternVertSTree(0,stree{chain2(iChain2),2},stree{p2,2},1/theta, fitLB, freqDistrCurr,eps,false);
        p = p2;
        addOrder(count) = chain2(iChain2);
        count = count + 1;
        iChain2 = iChain2 + 1;
    end
%     if curFit(p) > 0
        ub = ub + log(probMutEvent);
%     end
    if ub <= bestLikelihood
        ubCut = true;
        return;
    end
end

for i = iChain1:length(chain1)
    p = stree{chain1(i),3};
    fitLB(chain1(i)) =  fitLB(p) - freqCoeff(chain1(i))*theta/(n-(lastOrder+i-lastCh1)+1);
    [freqDistrCurr,probMutEvent,phi] = estimFreqInternVertSTree(0,stree{chain1(i),2},stree{p,2},1/theta, fitLB, freqDistrCurr,eps,false);
    addOrder(count) = chain1(i);
    count = count+1;
    %     if curFit(p) > 0
        ub = ub + log(probMutEvent);
%     end
%     ub = ub + log(probMutEvent);
    if ub <= bestLikelihood
        ubCut = true;
        return;
    end
end
for i = iChain2:length(chain2)
    p = stree{chain2(i),3};
    fitLB(chain2(i)) =  fitLB(p) - freqCoeff(chain2(i))*theta/(n-(lastOrder+i-lastCh2)+1);
    [freqDistrCurr,probMutEvent,phi] = estimFreqInternVertSTree(0,stree{chain2(i),2},stree{p,2},1/theta, fitLB, freqDistrCurr,eps,false);
    addOrder(count) = chain2(i);
    count = count + 1;
%     if curFit(p) > 0
        ub = ub + log(probMutEvent);
%     end
    if ub <= bestLikelihood
        ubCut = true;
        return;
    end
end