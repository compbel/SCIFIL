function [likelihood, fitInfer, orderMutInfer,stats] = fitBruteForceBackTrack2(stree,theta,fmax,fmin,freq,eps,Tmax,debugMode)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
degThr = 4;
root = find(cellfun(@isempty,stree(:,3)),1);
AM = streeOrder2AM(stree);
nHapl = size(AM,1);
G = digraph(AM);

if debugMode
    figure
    plot(G,'Layout','layered');
end

ord = flip(bfsearch(G,root))';
bif = zeros(1,nHapl);
nBif = 0;
for i = ord
    if length(stree{i,1}) >= 2
        nBif = nBif + 1;
        bif(nBif) = i;
    end
end
bif = bif(1:nBif);

for i = bif
    neigh_i = find(AM(i,:) > 0);
    
    chain_len = zeros(1,length(neigh_i));
    for j = 1:length(neigh_i)
        desc = bfsearch(G,neigh_i(j));
        chain_len(j) = length(desc);
    end
    hapl_neigh = cell2mat(stree(neigh_i,2));
    aux = [chain_len' freq(hapl_neigh) neigh_i'];
    aux = sortrows(aux,[1 2],{'descend' 'ascend'});
    neigh_i = aux(:,3);
    
    while length(neigh_i) >= 2
        desc1 = bfsearch(G,neigh_i(1));
        desc2 = bfsearch(G,neigh_i(2));
        desc = [i; desc1; desc2];
        nHapl_sub = length(desc);
        stree_sub = cell(2*length(desc),8);
        indDesc = zeros(1,nHapl);
        for j = 1:length(desc)
            indDesc(desc(j)) = j;
        end
        hapl_desc = cell2mat(stree(desc,2));
        freq_sub = freq(hapl_desc);
        freq_unnorm = freq_sub;
        if size(freq_sub,1) == 1
            freq_sub = freq_sub/sum(freq_sub,2);
        else
            freq_sub = freq_sub/sum(freq_sub,1);
        end
        for h = 1:length(desc)
            j = desc(h);
            stree_sub{h,2} = h;
            stree_sub{h,3} = indDesc(stree{j,3});
            if h == 1
                stree_sub{h,1} = [indDesc(neigh_i(1)) indDesc(neigh_i(2))];
            else
                stree_sub{h,1} = zeros(1,length(stree{j,1}));
                for k = 1:length(stree{j,1})
                    stree_sub{h,1}(k) = indDesc(stree{j,1}(k));
                end
            end
        end
        stree_sub{1,3} = [];
        for h = 1:length(desc)
            stree_sub{length(desc)+h,5} = freq_sub(h);
        end

%         AM_sub = streeOrder2AM(stree_sub);
%         G_sub = digraph(AM_sub);
%         figure
%         plot(G_sub,'Layout','layered');

    %     deg = sum(AM_sub(1,:),2);
    %     if deg >= degThr
    %         stree_sub = streeFreq2(stree_sub,freq_unnorm,degThr,eps);
    %         AM_sub = streeOrder2AM(stree_sub);
    %          G_sub = digraph(AM_sub);
    % %         figure
    % %         plot(G_sub,'Layout','layered');
    %     end

%         [likelihood_sub, fitInfer_sub, orderMutInfer_sub,stats_sub] = fitChains(stree_sub,nHapl_sub,theta,fmax,fmin,freq_sub,eps,Tmax);
        initOrder = [];
        initFit = [];
        initLikelihood = -Inf;
        [likelihood_sub, fitInfer_sub, orderMutInfer_sub,stats_sub] = fitChains_noRec(stree_sub,nHapl_sub,theta,fmax,fmin,freq_sub,eps,Tmax,initOrder,initFit,initLikelihood);
        if likelihood_sub == -Inf
            likelihood = -Inf;
            fitInfer = []; 
            orderMutInfer = [];
            stats = [];
            return
        end
        order_sub = orderMutInfer_sub;
%         order_sub = ones(1,length(desc));
%         try
%             order_sub(2:end) = orderMutInfer_sub + 1;
%         catch
%             ['ddd'];
%         end

        streeNew = stree;
        for j = 2:length(order_sub)
            u = desc(order_sub(j-1));
            v = desc(order_sub(j));
            streeNew{u,1} = v;
            streeNew{v,1} = [];
        end
        streeNew{i,1} = [streeNew{i,1} my_setdiff(stree{i,1},[neigh_i(1) neigh_i(2)])];

        stree = streeNew;
%         AM_new = streeOrder2AM(streeNew);
        AM_new = streeOrder2AM1(streeNew,nHapl);
        neigh_i = find(AM_new(i,:) > 0);
        
        G = digraph(AM_new);
%         for j = 1:length(neigh_i)
%              desc = bfsearch(G,neigh_i(j));
%              chain_len(j) = length(desc);
%         end
%         aux = [neigh_i' chain_len'];
%         aux = sortrows(aux,-2);
%         neigh_i = aux(:,1)';
%         if length(neigh_i) > 2
%             v = neigh_i(end);
%             neigh_i(end) = neigh_i(2);
%             neigh_i(2) = v;
%         end
        if debugMode
            figure
            plot(G,'Layout','layered');
        end
    end
end

% likelihood = likelihood_sub;
% fitInfer = fitInfer_sub; 
if ~isempty(bif)
    orderMutInfer = desc(order_sub);
    if length(orderMutInfer) < nHapl
        orderMutInfer = bfsearch(G,1);
    end
    stats = stats_sub;
else
    orderMutInfer = flip(ord)';
    stats = [];
end
orderMutInfer = orderMutInfer(2:end)-1;
[likelihood,fitInfer] = probTreeParamOrderSTreePartial(stree,orderMutInfer,theta,freq,eps,false);

