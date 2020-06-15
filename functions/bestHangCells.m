function M_best = bestHangCells(stree, M_obs, pFalseNeg, pFalsePos, pFalseNeg2, pFalsePos2)
    haplotypes = size(M_obs, 1);
    mutations = size(M_obs,2);
    M_best = zeros(haplotypes, mutations);
    tree_mut = (find(~cellfun(@isempty, stree(:,2)))).';
    prob_arr = zeros(1, length(tree_mut));
    for h=1:haplotypes
        best = -10e8;
        best_place = 1;
        
        for m=tree_mut
            prof = mutationProfile(stree,m,mutations);
            prob = probMatrHom(M_obs(h,:),prof,pFalseNeg,pFalsePos, pFalseNeg2, pFalsePos2);
            prob_arr(m) = prob;
            if (prob > best)
                best = prob;
                best_place = m;
            end
        end
        M_best(h,:) =  mutationProfile(stree,best_place,mutations);
    end
    k =1;
end