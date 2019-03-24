function lb = lowerBound1(stree,list, freq,curFit,nodes_count,eps,theta,fmax,fmin,currLikelihood)
% stree: [nextChild haplotype parent label frequency timet fitness oldChild]
lb = currLikelihood;
phi = curFit*freq';
% unvisited = setdiff(1:nodes_count,list);
unvisited = my_setdiff(1:nodes_count,list);
for v = unvisited
    p = stree{v+1,3};
    hapl_p = stree{p,2};
    if (curFit(hapl_p) > 0) && (curFit(hapl_p) <= phi)
        lb = lb - log(freq(hapl_p));
    end
end
