function [bestLikelihood, bestFit] = recEdgeLengthEstimationStree(step, stree, list, queue, curFit, theta, obsFreqLeafs,eps,Tmax, nodes_count, fmax, fmin, freqDistr, currLikelihood)
global gl_oter;
global l_min;
global bestOrder;
global filterFit;
global filterLB;
global filterPhi;
global filterFutureDeg;
bestLikelihood = -10e8;
bestFit = ones(1, length(obsFreqLeafs));
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
%     nIntern = sum(bintree(:,1) + bintree(:,2) > 0);
%     bestOrder = zeros(1,nIntern);
%     bestOrder = [];
if ~isempty(queue)
    for i = 1:length(queue)
        changed_fit = 0;
        newQueue = queue;
        j = length(list)+1;
        tv = (nodes_count-j+1)/theta;
        cr = 1;
        parent = stree{queue(i),3} ;
        %             % find cell child of a node
        %             for ch=stree{parent, 1}
        %                 if isempty(stree{ch,2})
        %                     cr = parent.children(ch).cell_index;
        %                     break;
        %                 end
        %             end
        try
            fr = stree{parent, 7};
        catch
            ['stop'];
        end
        Or = obsFreqLeafs(stree{parent,2});
        cl = queue(i);
        Ol = obsFreqLeafs(stree{cl,2});
        fl = fr - (log(eps/(1-eps))+log( Or) -log(Ol))/tv;
        if (fl > fmax || fl < fmin)
            filterFit = filterFit + 1;
            continue;
        end
        if isempty(curFit) % if we add first vertex to order after root we update fit vector with 2 found fitnesses
            changed_fit = 1;
            curFit = zeros(1, length(obsFreqLeafs));
            freqDistr = zeros(1, length(obsFreqLeafs));
            curFit(stree{1,2}) = 1;
            curFit(stree{queue(i),2}) = fl;
            freqDistr(stree{1,2}) = 1-eps;
            freqDistr(stree{queue(i),2}) = eps;
        end
        s_copy = stree;
        s_copy{cl,7} = fl;
        newFit = curFit;
        newFit(stree{cl,2}) = fl;
        s_copy{queue(i), 6} = tv;
        %             t_copy = tree.copy();
        %             t_copy.cells(cl).fit = fl;
        %             t_copy.nodes(queue(i)).time = tv;
        %             [list queue(i)]
        if ~isempty(list)
            %             currVert = find(freqDistr > 0); %???
            %                 for j =currVert
            %                     bt_copy(j, 7) = bt_copy(rightHapl(j),7);
            %                 end
            %                 bt_copy(queue(i),7) = bt_copy(rightHapl(queue(i)),7);
            %                 fit = zeros(1,length(tree.cells));
            %                 for c=1:length(tree.cells)
            %                     fit(c) = tree.cells(c).fit;
            %                 end
            [freqDistrNew,probMutEvent,phi] = estimFreqInternVertSTree(0, s_copy{queue(i), 2}, s_copy{parent, 2}, abs(s_copy{queue(i), 6}- s_copy{list(end)+1, 6}), newFit, freqDistr, eps,true);
            if currLikelihood - log(probMutEvent) > l_min
                filterLB = filterLB + 1;
                continue;
            end
%             if fl <= phi
%                 filterPhi = filterPhi + 1;
%                 continue;
%             end
            lb = lowerBound1(s_copy, [list queue(i)-1], freqDistrNew,newFit,nodes_count,eps,theta,fmax,fmin,currLikelihood - log(probMutEvent));
            if lb > l_min
                filterFutureDeg = filterFutureDeg + 1;
                continue;
            end
        else
            probMutEvent = 1;
            freqDistrNew = freqDistr;
%             [freqDistrNew,probMutEvent] = estimFreqInternVertSTree(0, s_copy{queue(i), 2}, s_copy{parent, 2}, abs(s_copy{queue(i), 6}- s_copy{parent, 6}), newFit, freqDistr, eps);
        end        
        
        newQueue(i) = [];
        %add all children mutations
        for c=s_copy{queue(i),1}
            newQueue = [newQueue c];
        end
%         newQueue1 = zeros(1,length(newQueue) + length(s_copy{queue(i),1})-1);
%         j = 1;
%         for ii = 1:length(newQueue)
%             if ii ~= i
%                 newQueue1(j) = newQueue(ii);
%                 j = j + 1;
%             end
%         end
%         for c=s_copy{queue(i),1}
%             newQueue1(j) = c;
%             j = j + 1;
%         end
%         newQueue = newQueue1;
        %             [list queue(i)]
        %             currLikelihood
        [l,b] = recEdgeLengthEstimationStree(step+1, s_copy, [list queue(i)-1], newQueue, newFit, theta, obsFreqLeafs,eps,Tmax, nodes_count, fmax, fmin, freqDistrNew, currLikelihood - log(probMutEvent));
        if -l < -bestLikelihood
            bestLikelihood = l;
            bestFit = b;
        end
        if changed_fit
            curFit = [];
            freqDistr = [];
        end
    end
else
    gl_oter = gl_oter +1;
%     [bestLikelihood,bestFit] = probTreeParamOrderSTree(stree,list,theta,obsFreqLeafs,eps);
    [bestLikelihood,bestFit] = probTreeParamOrderSTreePartial(stree,list,theta,obsFreqLeafs,eps,true);
    if -bestLikelihood < l_min
        l_min = -bestLikelihood;
        bestOrder = list;
        l_min
    end
    return
end