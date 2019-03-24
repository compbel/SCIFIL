function [l_rec, fitInfer1, orderMutInfer] = fitLocalSearch(stree,obsFreqLeafs,fmax,fmin,theta,eps,alpha,beta,M_noise) 
Tmax = 150000;
debugMode = false;
curTree = stree;
M_infer = stree2MutMatr(curTree);
nIntern = size(M_infer,2);

for iter=1:1
    [l_f_rec, fitInfer1, orderMutInfer,stats1] = fitBruteForceBackTrack2(curTree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax,debugMode);
    l_m_rec = probMatr(M_infer,M_noise,alpha, beta);
%     l_rec = l_f_rec + l_m_rec;
    l_rec = l_f_rec;
    
    improve = true;
    improveCount = 0;
    while improve
        improve = false;
        % parents(i,j) == 1 if j is a predecessor of i, including itself
        parents = zeros(nIntern+1,nIntern+1);
        for i=2:nIntern+1
            cur = i;
            while ~isempty(cur)
                parents(i,cur) = 1;
                cur = curTree{cur, 3};
            end
        end
        for f=2:nIntern+1
            for t=1:nIntern+1
                % don't move if from is a predecessor of to
                if parents(t,f) || curTree{f, 3} == t
                    continue;
                end
                new_stree = transformStree(curTree, f, t);
                %            AM = streeOrder2AM(new_stree);
                %            bg1 = biograph(AM,labels1);
                %             view(bg1)
                [l_f, fit_n, order_n,stats1] = fitBruteForceBackTrack2(new_stree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax,debugMode);
                M_infer = stree2MutMatr(new_stree);
                l_m = probMatr(M_infer,M_noise,alpha, beta);
%                 l_n = l_f + l_m;
                l_n = l_f;
                if (l_n > l_rec)
                    l_n;
                    l_f_rec = l_f;
                    l_m_rec = l_m;
                    l_rec = l_n;
                    best_fit_iter = fit_n;
                    best_order_iter = order_n;
                    best_tree_iter = new_stree;
                    improve = true;
                end
            end
        end
        if improve
            curTree = best_tree_iter;
            fitInfer1 = best_fit_iter;
            orderMutInfer = best_order_iter;
            improveCount = improveCount + 1;
            if improveCount >= 1
                return;
            end
        end
    end
    l_rec
    
end



