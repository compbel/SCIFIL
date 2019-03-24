function AM = streeOrder2AM(stree)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
intern = find(~cellfun(@isempty,stree(:,2)));
nIntern = length(intern);
AM = zeros(nIntern,nIntern);
for i = intern'
    for j = intern'
        if ismember(j,stree{i,1})
            AM(i,j) = 1;
        end
    end
end