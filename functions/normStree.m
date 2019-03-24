function [newStree, newFreq, indLeft] = normStree(stree,freq,eps)
% stree: [nextChild haplotype parent label frequency timet fitness oldChildren]
isLeft = zeros(1,length(freq));
isNoCell = zeros(1,length(freq));
newFreq = zeros(size(freq,1),size(freq,2));
newStree = stree;
for i = 1:length(freq)
    if isempty(stree{i,2})
        isNoCell(i) = 1;
    else
       freq_node = freq(stree{i,2});
       [m,j] = max(freq_node);
       isLeft(stree{i,2}(j)) = 1;
       newFreq(stree{i,2}(j)) = sum(freq_node,1);
       newStree{i,2} = stree{i,2}(j);
    end
end
empMut = find(isNoCell == 1);
empHapl = find(isLeft == 0);

for i = 1:length(empMut)
    newStree{empMut(i),2} = empHapl(i);
    newFreq(empHapl(i)) = eps;
end

indLeft = find(isLeft == 1);
