function stree = AMInfScite2STree(AM)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
stree = cell(2*length(AM),7);
for i = 1:length(AM)
    p = find(AM(:,i) > 0);
    if ~isempty(p)
        stree{i,3} = p;
    end
    stree{length(AM)+i,3} = i;
    ch = find(AM(i,:) > 0);
    stree{i,1} = ch;
    stree{i,2} = i;
end