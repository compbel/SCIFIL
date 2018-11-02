function bintree = estimTime1(bintree,Tmax,eps) 
% bintree: [leftChild righChild parent label frequency timet fitness]
zerofreq = (bintree(:,5) == 0);
bintree(zerofreq,5) = eps^2;
internal = (find(bintree(:,1)+bintree(:,2) > 0))';
for j = length(internal):-1:1
    p = internal(j);
    f2 = bintree(p,7);
    f1 = bintree(bintree(p,1),7);
    o2_all = bintree(p,5);
    o1_all = bintree(bintree(p,1),5);
    o1 = o1_all/(o1_all+o2_all);
    o2 = o2_all/(o1_all+o2_all);
    Tmin = max(bintree(bintree(p,1),6),bintree(bintree(p,2),6))+1;
%         Tmrca = getTMRCA(f1,o1,f2,o2,Tmax,eps);
%         if Tmin > Tmax
%             Tmax = 2*Tmax;
%         end
    Tmrca = getTMRCA2(f1,o1,f2,o2,Tmin,Tmax,eps);
    bintree(p,6) = Tmrca;
end

% AM = treeAMdirect(bintree);
% nBinNodes = size(bintree,1);
% labels = cell(1,nBinNodes);
% for i = 1:nBinNodes
%     labels{i} = ['v' int2str(i) '_l=' int2str(bintree(i,4)) '_t=' num2str(bintree(i,6))];
% end
% bg = biograph(AM,labels);
% view(bg)

T = max(bintree(:,6));
bintree(:,6) = T - bintree(:,6);




