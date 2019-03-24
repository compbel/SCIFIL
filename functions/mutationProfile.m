function mut = mutationProfile(stree, mutation, m)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
    mut = zeros(1,m);
    cur = mutation;
    while (~isempty(stree{cur, 3}))
        mut(stree{cur, 2}-1) = 1;
        cur = stree{cur, 3};
    end
end