function AM = streeOrder2AM1(stree,nIntern)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
% intern = find(~cellfun(@isempty,stree(:,2)));
% nIntern = length(intern);
AM = zeros(nIntern,nIntern);
for i = 1:nIntern
    AM(i,stree{i,1}) = 1;
end