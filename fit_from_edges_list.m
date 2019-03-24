disp('Starting SCIFIL')
addpath('functions');
if (~exist('method', 'var'))
    method = 'brute_force';
end

if (~exist('output', 'var'))
    output = 'out.txt';
end



M_noise = dlmread(matrix_file); % matrix_file = 'data/scite/dataHou18.csv'

M_noise = M_noise';

n = size(M_noise,1);
m = size(M_noise,2);
if (~exist('nRep', 'var'))
    nRep = 0;
end
m = m+nRep;

AMscite = scite2Tree(gv_file,n,m); % gv_file = 'data/scite/dataHou18_map0.gv'
stree = AMInfScite2STree(AMscite(1:(m+1),1:(m+1)));

cellCounts = zeros(m+1,1);
for i = 1:n
    ind = find(AMscite(:,m+1+i) > 0);
    cellCounts(ind) = cellCounts(ind) + 1/length(ind);
    %         cellCounts(ind) = cellCounts(ind) + 1;
end
obsFreqLeafs = cellCounts/sum(cellCounts,1);

% [stree, obsFreqLeafs, indLeft] = normStree(stree,obsFreqLeafs,eps1);

fmax = 1.5;
fmin = 1;
theta = 0.01;
eps = 1e-5;
Tmax = 150000;
timeLimit = 600;
alpha = 0.2;
beta = 0.00001;
debugMode = false;


% labels1 = cell(1,n+m+1);
% for i = 1:(m+1)
%     labels1{i} = ['m',int2str(i-1)];
% end
% for i =1:n
%     labels1{m+1+i} = ['c',int2str(i),'_'];%,num2str(obsFreqLeafs(i))];
% end
% bg1 = biograph(AMscite,labels1);
% view(bg1)

% intern = find(bintree(:,1) + bintree(:,2) > 0);
% orderTrue = intern;
% [leftHapl,rightHapl] = getProducedHapl(bintree);
% [likelihood2,fitInfer1] = probTreeParamOrder(bintree,orderTrue,leftHapl,rightHapl,theta,obsFreqLeafs,eps);

% [likelihood,fit,order,stats] = bintreeFitBruteForce(bintree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax,timeLimit);
% [likelihood, fit, order] = fitLocalSearch(stree,obsFreqLeafs,5*fmax,0.5*fmin,theta,eps,alpha,beta,M_noise)
% fill all fits according to initial order of mutations

if (strcmp(method, 'brute_force'))
    disp('Starting heuristic estimation');
    [likelihood,fit, order,stats] = fitBruteForceStree(stree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax);
else
    disp('Starting brute force calculation');
    [likelihood, fit, order,stats] = fitBruteForceBackTrack2(stree,theta,5*fmax,0.5*fmin,obsFreqLeafs,eps,Tmax,debugMode);
end
true_fit_order = fit;
disp("")
disp("Found fitness vector:");
disp(true_fit_order);
disp(['writeing results to ' output])
fileID = fopen(output,'w');
fprintf(fileID,'%6.4f ', true_fit_order);
fprintf(fileID,'\n');
fclose(fileID);
