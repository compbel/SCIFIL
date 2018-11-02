function stree = AMtoSTree(AM, nIntern, nUnCells)
stree = cell(nIntern+nUnCells+1,7);
for i = 1:(nIntern+nUnCells+1)
    p = find(AM(:,i) > 0);
    if ~isempty(p)
        stree{i,3} = p;
    end
    ch = find(AM(i,1:(nIntern+1)) > 0);
    stree{i,1} = ch;
    if i <= nIntern + 1
        h = find(AM(i,(nIntern+2):end) > 0);
        if ~isempty(h)
            stree{i,2} = h;
        end
        stree{i,4} = i-1;
    end
end