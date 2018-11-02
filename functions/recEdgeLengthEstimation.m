function [bestLikelihood, bestFit] = recEdgeLengthEstimation(step, bintree, list, queue, theta, obsFreqLeafs,eps,Tmax, leftHapl, rightHapl, nodes_count, fmax, fmin, freqDistr, currLikelihood)  
    global gl_oter;   
    global l_min;
    global bestOrder;
    global filterFit;
    global filterLB;
    global time_limit;
    bestLikelihood = -10e8;
    bestFit = ones(1, length(obsFreqLeafs));
    c = toc;
    if time_limit < c
        return
    end    
%     nIntern = sum(bintree(:,1) + bintree(:,2) > 0);
%     bestOrder = zeros(1,nIntern);
%     bestOrder = [];
    if ~isempty(queue)
        for i = 1:length(queue)
            newQueue = queue;            
            j = length(list)+1;
            tv = (nodes_count-j+1)/theta;
            fr = bintree(rightHapl(queue(i)),7);
            Or = bintree(rightHapl(queue(i)), 5);
            Ol = bintree(leftHapl(queue(i)), 5);
            fl = fr - (log(eps/(1-eps))+log( Or) -log(Ol))/tv;
            if (fl > fmax || fl < fmin)
                filterFit = filterFit + 1;
                continue;
            end
            bt_copy = bintree;
            bt_copy(leftHapl(queue(i)), 7) = fl;
            bt_copy(queue(i), 6) = tv;
%             [list queue(i)]
            if ~isempty(list)
                currVert = find(freqDistr > 0);
                for j =currVert
                    bt_copy(j, 7) = bt_copy(rightHapl(j),7);
                end
                bt_copy(queue(i),7) = bt_copy(rightHapl(queue(i)),7);
                [freqDistrNew,probMutEvent] = estimFreqInternVert(bt_copy,queue(i),list(end),freqDistr,eps);
                if currLikelihood - log(probMutEvent) > l_min
                    filterLB = filterLB + 1;
                    continue;
                end
            else
                probMutEvent = 1;
                freqDistrNew = freqDistr;
            end
           

            newQueue(i) = [];
            % if left child is not leaf, than add it to queue
            if bintree(bintree(queue(i), 1), 1) ~= 0
                newQueue = [newQueue bintree(queue(i), 1)];
            end
            % if right child is not leaf, than add it to queue
            if bintree(bintree(queue(i), 2), 1) ~= 0
                newQueue = [newQueue bintree(queue(i), 2)];
            end
%             [list queue(i)]
%             currLikelihood
            [l,b] = recEdgeLengthEstimation(step+1, bt_copy, [list queue(i)], newQueue, theta, obsFreqLeafs,eps,Tmax, leftHapl, rightHapl, nodes_count, fmax, fmin, freqDistrNew, currLikelihood - log(probMutEvent));
            if -l < -bestLikelihood
                bestLikelihood = l;
                bestFit = b;
            end
        end
    else
        gl_oter = gl_oter +1;
        [bestLikelihood,bestFit] = probTreeParamOrder(bintree,list,leftHapl,rightHapl,theta,obsFreqLeafs,eps);
        if -bestLikelihood < l_min
            l_min = -bestLikelihood;
            bestOrder = list;
            l_min
        end
        return 
    end