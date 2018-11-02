function streeInfer = estimTimeSTree(stree,fit,obsFreqLeafs,internOrder,Tmax,eps) 
% stree: [nextChild haplotype parent label frequency timet fitness] row 1
% is mutation 0

nIntern = length(internOrder);
timesInfer = zeros(1,nIntern);
for i = length(internOrder):-1:1
    j = internOrder(i) +1; % +1???
    p = stree{j,3};
    f1 = fit(stree{j,2});
    f2 = fit(stree{p,2});
    o2 = obsFreqLeafs(stree{p,2});
    o1 = obsFreqLeafs(stree{j,2});
    if ~isempty(stree{j,1})
        Tmin = max(timesInfer(stree{j,1}-1))+1;
    else
        Tmin = 1;
    end
%         Tmrca = getTMRCA(f1,o1,f2,o2,Tmax,eps);
%         if Tmin > Tmax
%             Tmax = 2*Tmax;
%         end
    Tmrca = getTMRCA2(f1,o1,f2,o2,Tmin,Tmax,eps);
    timesInfer(j-1) = Tmrca;
end

% AM = treeAMdirect(bintree);
% nBinNodes = size(bintree,1);
% labels = cell(1,nBinNodes);
% for i = 1:nBinNodes
%     labels{i} = ['v' int2str(i) '_l=' int2str(bintree(i,4)) '_t=' num2str(bintree(i,6))];
% end
% bg = biograph(AM,labels);
% view(bg)

T = max(timesInfer);
timesInfer = T - timesInfer;
streeInfer = stree;
for i = 1:nIntern
    streeInfer{i+1,6} = timesInfer(i);
end
for i = (nIntern+2):size(stree,1)
    streeInfer{i,6} = T;
end




