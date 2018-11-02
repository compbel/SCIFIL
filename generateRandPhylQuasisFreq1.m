function [M,AM,cellDistr,timesMut,t1,bintree,stree,popDynam,order,orderMut] = generateRandPhylQuasisFreq1(m,theta,f,eps,randTimeInt,randMut)
% tree: [leftChild righChild parent label cellCount]
% bintree: [leftChild righChild parent label frequency timet fitness]
% stree: [nextChild haplotype parent label frequency timet fitness]
% root - 1st row; edges labelled by their second ends

bintree = [0 0 0 0 1 0 1];
leafs = [1];
nNodes = 1;
nodeMuts = cell(1,1);
nodeMuts{1} = [];
mutAppear = zeros(1,m);
t = 0;
dt = 0.001;
expTime = 1/theta;

leaf = leafs(1);
mut = randsample(setdiff(1:m,nodeMuts{leaf}),1);
mutAppear(mut) = 1;
bintree(leaf,1) = nNodes + 1;
bintree(leaf,2) = nNodes + 2;
bintree(leaf,4) = mut;
bintree(leaf,6) = t;
fit = 1 + rand*(f-1);
bintree = [bintree; [0 0 leaf 0 eps 0 fit]; [0 0 leaf 0 bintree(leaf,5)-eps 0 bintree(leaf,7)]];
nodeMuts{nNodes + 1} = [nodeMuts{leaf} mut];
nodeMuts{nNodes + 2} = nodeMuts{leaf};
leafs(leafs == leaf) = [];
leafs = [leafs nNodes+1 nNodes+2];
nNodes = nNodes + 2;
tSincePrevMut = 0;
if randTimeInt == true
    Time = round(poissrnd(expTime)); 
else
    Time = expTime;
end

popDynam = [];
% popDynam = zeros(length(0:dt:(expTime*m)),m+1);
% ipop = 1;

while sum(mutAppear,2) < m
% for t = 1:T
%     theta = N*m*mutRate;
%     prob = exprnd(theta);
    p = rand;
%     if p < theta
    if tSincePrevMut == Time
%         freqDistr = zeros(1,length(leafs));
%         for j = 1:length(leafs)
%             l = leafs(j);
%             freqDistr(j) = bintree(l,5)*exp(bintree(l,7)*tSincePrevMut);
%         end
%         freqDistr = freqDistr/sum(freqDistr,2);
        
       xCurr = bintree(leafs,5);
%        popDynam(ipop,leafs) = xCurr;
       for tau = dt:dt:tSincePrevMut
           phi = xCurr'*bintree(leafs,7);
           xNext = xCurr + dt*(bintree(leafs,7).*xCurr - phi*xCurr);
           xCurr = xNext;
%            ipop = ipop + 1;
%            popDynam(ipop,leafs) = xCurr;
       end
       freqDistr = xCurr';
       bintree(leafs,5) = freqDistr';
%        phi = xCurr'*bintree(leafs,7);
        
        if randMut == true
            CDF = [0 cumsum(freqDistr)];
            p1 = rand;
            j = find(CDF >= p1,1)-1;
        else
            [mf,j] = max(freqDistr);
        end
        leaf = leafs(j);
        
        possMuts = setdiff(1:m,union(nodeMuts{leaf},find(mutAppear)));
        mut = possMuts(randi(length(possMuts)));
%         mut = randsample(possMuts,1);
%         if mutAppear(mut) == 1
%             p_rep = rand;
%             if p_rep < probRep
%                 continue;
%             end
%         end
        mutAppear(mut) = 1;
        bintree(leaf,1) = nNodes + 1;
        bintree(leaf,2) = nNodes + 2;
        bintree(leaf,4) = mut;
        bintree(leaf,6) = t;
        
        fit = 1 + rand*(f-1);
%         fit = phi + rand*(f-phi);
%         bintree = [bintree; [0 0 leaf 0 eps 0 fit]; [0 0 leaf 0 bintree(leaf,5)-eps 0 bintree(leaf,7)]];
        bintree = [bintree; [0 0 leaf 0 eps*bintree(leaf,5) 0 fit]; [0 0 leaf 0 (1-eps)*bintree(leaf,5) 0 bintree(leaf,7)]];

        nodeMuts{nNodes + 1} = [nodeMuts{leaf} mut];
        nodeMuts{nNodes + 2} = nodeMuts{leaf};
        leafs(leafs == leaf) = [];
        leafs = [leafs nNodes+1 nNodes+2];
        nNodes = nNodes + 2;
        tSincePrevMut = 0;
        if randTimeInt == true
            Time = round(poissrnd(expTime)); 
        else
            Time = expTime;
        end
    else
        tSincePrevMut = tSincePrevMut + 1;
        t = t+1
%         for j = leafs
%             bintree(j,5) = bintree(j,7)*bintree(j,5);
%         end
%         bintree(leafs,5) = bintree(leafs,5)/sum(bintree(leafs,5),1);

%         phi = sum(bintree(leafs,7).*bintree(leafs,5),1);
%         for j = leafs
%             bintree(j,5) = bintree(j,5) + bintree(j,7)*bintree(j,5) - phi*bintree(j,5);
%          end
    end
end

intern = setdiff((1:size(bintree,1)),leafs);
nIntern = size(intern,2);
% T_after = round((max(bintree(intern,6))-min(bintree(intern,6)))/(size(intern,2))-1);
T_after = expTime;

% for j = 1:length(leafs)
%     l = leafs(j);
%     freqDistr(j) = bintree(l,5)*exp(bintree(l,7)*tSincePrevMut);
% end
% freqDistr = freqDistr/sum(freqDistr,2);

xCurr = bintree(leafs,5);
% popDynam(ipop,leafs) = xCurr;
for tau = dt:dt:T_after
    phi = xCurr'*bintree(leafs,7);
    xNext = xCurr + dt*(bintree(leafs,7).*xCurr - phi*xCurr);
    xCurr = xNext;
%     ipop = ipop + 1;
%     popDynam(ipop,leafs) = xCurr;
end
freqDistr = xCurr';
bintree(leafs,5) = freqDistr';
       
t1 = t+T_after;

freqDistr = bintree(leafs,5);
cellDistr = freqDistr;
nUnCells = size(leafs,2);
M = zeros(nUnCells,nIntern);
c=1;
for i = 1:size(cellDistr,1)
    if cellDistr(i) > 0
        for k = nodeMuts{leafs(i)}
            M(c,k) = 1;
        end
        bintree(leafs(i),4) = nIntern + c;
        c = c+1;
    end
end

% AM_bin = treeAMdirect(bintree);
% labels = cell(1,nNodes);
% for i = 1:nNodes
%     labels{i} = ['v' int2str(i) '_' int2str(bintree(i,4))];
% end
% bg = biograph(AM_bin,labels);
% view(bg)

timesMut = zeros(1,nIntern);
AM = zeros(nIntern+nUnCells+1,nIntern+nUnCells+1);
for i = 1:size(bintree,1)
    if ismember(i,leafs)&& (cellDistr(leafs==i)== 0)
        continue;
    end
    v = i;
    par = bintree(v,3);
    if ~ismember(i,leafs)
        timesMut(bintree(i,4)) = bintree(i,6);
    end
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

% stree: [nextChild haplotype parent label frequency timet fitness]
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

% labels1 = cell(1,nIntern+nUnCells+1);
% for i = 1:(nIntern+1)
%     labels1{i} = ['m',int2str(i-1)];
% end
% for i =1:nUnCells
%     labels1{nIntern+1+i} = ['c',int2str(i),'_',int2str(cellCounts(i))];
% end
% % bg1 = biograph(AM(1:(nIntern+1),(1:nIntern+1)),labels1(1:(nIntern+1)));
% bg1 = biograph(AM,labels1);
% view(bg1)
aux = [(1:nIntern)' bintree(intern,6)];
aux = sortrows(aux,2);
order = intern(aux(:,1));
orderMut = bintree(order,4);
['Done']