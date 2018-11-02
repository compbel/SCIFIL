function stree = bintreeToStree(bintree, nIntern,nUnCells, cellDistr)
AM = zeros(nIntern+nUnCells+1,nIntern+nUnCells+1);
leafs = (find(bintree(:,1)+bintree(:,2) == 0))';
for i = 1:size(bintree,1)
    if ismember(i,leafs)&& (cellDistr(leafs==i)== 0)
        continue;
    end
    v = i;
    par = bintree(v,3);
    while (par~=0)&&(bintree(par,1)~=v)
        v = par;
        par = bintree(v,3);
    end
    if par > 0
        AM(bintree(par,4)+1,bintree(i,4)+1) = 1;
    else
        AM(1,bintree(i,4)+1) = 1;
    end
end

labels1 = cell(1,nIntern+nUnCells+1);
for i = 1:(nIntern+1)
    labels1{i} = ['m',int2str(i-1)];
end
for i =1:nUnCells
    labels1{nIntern+1+i} = ['c',int2str(i),'_',num2str(cellDistr(i))];
end
bg1 = biograph(AM,labels1);
view(bg1)
stree = AMtoSTree(AM,nIntern,nUnCells);