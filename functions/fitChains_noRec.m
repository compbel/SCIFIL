function [likelihood, fitInfer, orderMutInfer,stats] = fitChains_noRec(stree,nHapl,theta,fmax,fmin,obsFreqLeafs,eps,Tmax,initOrder,initFit,initLikelihood)

gl_oter = 0;
bestOrder = initOrder;
bestFit = initFit;
bestLikelihood = initLikelihood;
filterFit = 0;
filterLB = 0;
filterPhi = 0;
avCutPoint = 0;
filterFutureDeg = 0;
times = zeros(1,3);

root = find(cellfun(@isempty,stree(:,3)),1);
AM = streeOrder2AM1(stree,nHapl);
G = digraph(AM);
chains = cell(1,2);
children = stree{root,1};
for i = 1:length(children)
    chains{i} = bfsearch(G,children(i));
end

% figure
% plot(G,'Layout','layered');

n = nHapl-1;
t_last = n/theta;
r = min(length(chains{1}),length(chains{2}));
if length(chains{1}) > r
    ch = chains{1};
    chains{1} = chains{2};
    chains{2} = ch;
end
subset = zeros(1,r);
freqDistr = zeros(1,nHapl);
freqDistr(root) = 1;
currFit = zeros(1,nHapl);
currFit(root) = 1;
stree{root,6} = 0;

freqCoeff = zeros(1,nHapl);
for i = 2:nHapl
    parent = stree{i,3} ;
    Op = obsFreqLeafs(stree{parent,2});
    Ov = obsFreqLeafs(stree{i,2});
    freqCoeff(i) = log(eps/(1-eps))+log(Op)-log(Ov);
end

chain1 = chains{1}';
chain2 = chains{2}';
chainsCharVect = zeros(1,nHapl);
chainsCharVect(chain1) = 1;
chainsCharVect(chain2) = 2;

stackElement = zeros(1,n*(n+1)/2);
stackPosition = zeros(1,n*(n+1)/2);
stackElement(1:(n-r+1)) = (n-r+1):-1:1;
stackPosition(1:(n-r+1)) = 1;

stackCurrOrder = cell(1,n*(n+1)/2);
stackCurrFit = cell(1,n*(n+1)/2);
stackCurrLikelihood = zeros(1,n*(n+1)/2);
stackCurrElem2 = zeros(1,n*(n+1)/2);
stackFreqDistr = cell(1,n*(n+1)/2);
% stackMaxFit = zeros(1,n*(n+1)/2);


for i = 1:(n-r+1)
    stackCurrFit{i} = currFit;
    stackFreqDistr{i} = freqDistr;
end

top = n-r+1;
while top  > 0
    elem = stackElement(top);
    pos = stackPosition(top);
    currOrder = stackCurrOrder{top};
    currElem2 = stackCurrElem2(top);
    currLikelihood = stackCurrLikelihood(top);
    currFit = stackCurrFit{top};
    freqDistr = stackFreqDistr{top};
%     fm = stackMaxFit(top);
    top = top - 1;
    subset(pos) = elem;
 
    if pos >=2
        gap = subset(pos) - subset(pos-1) - 1;
    else
        gap = subset(pos)-1;
    end

%     orderAppend = [chain2((currElem2+1):(currElem2+gap)) chain1(pos)];
    newOrder = [currOrder chain2((currElem2+1):(currElem2+gap)) chain1(pos)];
    currElem2new = currElem2+gap;
    if subset(end) ~= 0
        gap1 = n - subset(end);
%         orderAppend = [orderAppend chain2((currElem2new+1):(currElem2new+gap1))];
        newOrder = [newOrder chain2((currElem2new+1):(currElem2new+gap1))];
    end
%     newOrder = [currOrder orderAppend];
    newLikelihood = currLikelihood;
    % time = toc;
    % times(2) = times(2) + time;
    % 
    % tic
    newFit = currFit;
    freqDistrPrev = freqDistr;
    filterPassed = true;
    
        
%     if (subset(3) == 15) && (subset(4) == 17)
%        ['stop']
%     end

%     if (subset(1) == 2)
%        ['stop']
%     end
    
    
%     for i = (length(currOrder)+1):(length(currOrder) + length(orderAppend))
    for i = (length(currOrder)+1):(length(newOrder))
        v = newOrder(i);
        tv = (i-1)/theta;
        parent = stree{v,3} ;
        fp = newFit(parent);
%         Op = obsFreqLeafs(stree{parent,2});
%         Ov = obsFreqLeafs(stree{v,2});
%         fv = fp - (log(eps/(1-eps))+log(Op)-log(Ov))/(t_last-tv);
        fv = fp - freqCoeff(v)/(t_last-tv);
        newFit(stree{v,2}) = fv;
        stree{v,6} = tv;
        if i == 1
            timeInterv = 0;
        else
            prev = newOrder(i-1);
            timeInterv = stree{v,6}-stree{prev,6};
        end
        [freqDistrNew,probMutEvent,phi] = estimFreqInternVertSTree(0,stree{v,2},stree{parent,2},timeInterv, newFit, freqDistrPrev, eps,false);
%         [freqDistrNew,probMutEvent] = estimFreqInternVertSTree_norec(0,stree{v,2},stree{parent,2},timeInterv, newFit, freqDistrPrev, eps);

        freqDistrPrev = freqDistrNew;
        newLikelihood = newLikelihood + log(probMutEvent);
        if newLikelihood < bestLikelihood
            filterLB = filterLB + 1;
            filterPassed = false;
            avCutPoint = avCutPoint + i;
            break;
        end
%         if fv <= phi
%             filterPhi = filterPhi + 1;
%             filterPassed = false;
%             break;
%         end
%         if (fv > fmax || fv < fmin)
%             filterFit = filterFit + 1;
%             filterPassed = false;
%             break;
%         end
%         ub = upperBound1(stree, newOrder(1:i), freqDistrNew,newFit,n,eps,theta,fmax,fmin,newLikelihood);
%         if ub < bestLikelihood
%             filterFutureDeg = filterFutureDeg + 1;
%             filterPassed = false;
%             break;
%         end
%         if i == length(newOrder)
%             lastCh1 = pos;
%         else
%             lastCh1 = pos - 1;
%         end
%         if chainsCharVect(v) == 1
%             lastCh1 = pos;
%         else
%             lastCh1 = pos - 1;
%         end
%         [ubCut,ub] = upperBound3(stree,newFit,newLikelihood,chain1,chain2,lastCh1,i-lastCh1,i,freqCoeff,theta,freqDistrNew,bestLikelihood,eps);
%         if ubCut
%             filterLB = filterLB + 1;
%             filterPassed = false;
%             avCutPoint = avCutPoint + i;
%             break;
%         end
    end
    
    if pos < r
        [ubCut,ub,addOrder,fitLB] = upperBound3(stree,newFit,newLikelihood,chain1,chain2,pos,length(newOrder)-pos,i,freqCoeff,theta,freqDistrNew,bestLikelihood,eps);
        if ubCut
            filterLB = filterLB + 1;
            filterPassed = false;
            avCutPoint = avCutPoint + length(newOrder);
        else
            bestLikelihood = ub;
            bestFit = fitLB;
            bestOrder = [newOrder addOrder];
        end
        if filterPassed 
            for i = (n-r+pos+1):-1:(elem+1)
                top = top + 1;
                stackElement(top) = i;
                stackPosition(top) = pos+1;
                stackCurrOrder{top} = newOrder;
                stackCurrFit{top} = newFit;
                stackCurrLikelihood(top) = newLikelihood;
                stackCurrElem2(top) = currElem2new;
                stackFreqDistr{top} = freqDistrPrev;
%                 stackMaxFit(top) = fm_new;
            end
        end
    else
        if (newLikelihood > bestLikelihood) && (filterPassed)
            bestLikelihood = newLikelihood;
            bestFit = newFit;
            bestOrder = newOrder;
            gl_oter = gl_oter + 1;
       end
    end
end

likelihood = bestLikelihood;
fitInfer = bestFit;
orderMutInfer = [root bestOrder];
stats = [gl_oter filterFit filterLB filterPhi filterFutureDeg avCutPoint/filterLB times(1) times(2) times(3)];
