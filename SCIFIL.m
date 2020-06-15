% gv_file='data/dataHou18_map0_rep3.gv';
% names_file='data/dataHou18names.txt';
% nRep=3;
% n=58;
% m=18;
disp('Starting SCIFIL')
addpath('functions');
if (~exist('method', 'var'))
    method = 'heuristic';
end

if (~exist('output', 'var'))
    output = 'out.txt';
end

% M_noise = dlmread(matrix_file); % matrix_file = 'data/dataHou18.csv'
%
% M_noise = M_noise';
% 
% n = size(M_noise,1);
% m = size(M_noise,2);

% number of repeated mutation


if (exist('names_file', 'var'))
    fid = fopen( names_file); % 'data/dataHou18names.txt'
    names = textscan(fid,'%s');
    names = names{1,1};
    names = [{'empty'}; names];
else
    names = num2cell(1:m);
    names=cellstr(cellfun(@num2str,names,'uni',false));
    names = names';
    names = [{'empty'}; names];
end

if (~exist('nRep', 'var'))
    nRep = 0;
end
if (~exist('noRep', 'var'))
    noRep = 0;
end
if nRep > 0
    names{end+1} = [names{nRep+1} '_1'];
    noRep = 1;
end
m = m+noRep;



AMscite = scite2Tree(gv_file,n,m); % gv_file = 'data/dataHou18_map0_rep3.gv'
stree = AMInfScite2STree(AMscite(1:(m+1),1:(m+1)));

cellCounts = zeros(m+1,1);
for i = 1:n
    ind = find(AMscite(:,m+1+i) > 0);
    cellCounts(ind) = cellCounts(ind) + 1/length(ind);
    %         cellCounts(ind) = cellCounts(ind) + 1;
end
cellCounts(cellCounts==0) = 0.1;
obsFreqLeafs = cellCounts/sum(cellCounts,1);

% [stree, obsFreqLeafs, indLeft] = normStree(stree,obsFreqLeafs,eps1);

fmax = 1.5;
fmin = 1;
if (~exist('theta', 'var'))
    theta = 0.01;
end
eps = 1e-4;
Tmax = 150000;
timeLimit = 600;
alpha = 0.2;
alpha2 = alpha;
beta = 0.00001;
beta2 = alpha*beta/2;
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
% [likelihood, fit, order] = fitLocalSearch(stree,obsFreqLeafs,5*fmax,0.5*fmin,theta,eps,alpha,beta,alpha2, beta2, M_noise)
% fill all fits according to initial order of mutations

if (strcmp(method, 'brute_force'))
    disp('Starting heuristic estimation');
    [likelihood,fitInfer, orderInfer,stats] = fitBruteForceStree(stree,theta,fmax,fmin,obsFreqLeafs,eps,Tmax);
else
    disp('Starting brute force calculation');
    [likelihood, fitInfer, orderInfer,stats] = fitBruteForceBackTrack2(stree,theta,5*fmax,0.5*fmin,obsFreqLeafs,eps,Tmax,debugMode);
end
disp("")
disp("Found fitness vector:");
disp(fitInfer);
disp(['writeing results to ' output])
fileID = fopen(output,'w');
fprintf(fileID,'%6.4f ', fitInfer);
fprintf(fileID,'\n');
fclose(fileID);


%         plot fitness landscape

labelsScite = cell(1,m+1+n);
for i = 0:m
    %         labelsScite{i+1} = ['m' int2str(i+1)];
    labelsScite{i+1} = names{i+1};
end
for i = 1:n
    labelsScite{m+1+i} = ['c' int2str(i)];
end
Gscite = digraph(AMscite(1:(m+1),1:(m+1)));
figure('Name','Mutation tree')
plot(Gscite,'NodeLabel',labelsScite(1:(m+1)));

AM_stree = AMscite(1:(m+1),1:(m+1));

AM = (AM_stree + AM_stree' > 0);
G = graph(AM);
DM_tree = distances(G,'Method','unweighted');
nHapl = size(DM_tree,1);

order = [1; orderInfer + 1];
appDist = zeros(nHapl,nHapl);
for i = 1:nHapl
    for j = 1:nHapl
        o1 = find(order == i);
        o2 = find(order == j);
        appDist(i,j) = abs(o1 - o2);
    end
end

fitNodes = zeros(1,nHapl);
for i = 1:nHapl
    fitNodes(i) = fitInfer(stree{i,2});
end

DM=DM_tree + 0.5*appDist;
DataMDS=cmdscale(DM,2);

maxCoord = max(max(abs(DataMDS)));
coeff = 5;
coeff1 = 0.1;
DataMDS = coeff*DataMDS/maxCoord;
lim1x = min(DataMDS(:,1)) -coeff1;
lim2x = max(DataMDS(:,1)) + coeff1;
lim1y = min(DataMDS(:,2)) -coeff1;
lim2y = max(DataMDS(:,2)) + coeff1;
[xq,yq] = meshgrid(lim1x:0.01:lim2x, lim1y:0.01:lim2y);
vq = griddata(DataMDS(:,1),DataMDS(:,2),fitNodes',xq,yq,'v4');
color = [0 0 0];
colors = ones(nHapl,3);

figure('Name','Fitness landscape')

contourf(xq,yq,vq,20)
axis off;
colormap jet
caxis([0.9 max(fitNodes)+0.3])
% caxis([0.9 1.4])
hold on
scatter(DataMDS(:,1),DataMDS(:,2),50,colors,'filled','MarkerEdgeColor','black')
for i=1:nHapl
    for j = 1:nHapl
        if AM_stree(i,j) == 1
            arrow(DataMDS(i,:), DataMDS(j,:),'Length', 8);
        end
        if AM_stree(j,i) == 1
            arrow(DataMDS(j,:), DataMDS(i,:),'Length', 8);
        end
    end
end
xoff = 0.1;
yoff= 0.05;
text(DataMDS(:,1)+xoff,DataMDS(:,2)+yoff, names, 'Fontsize', 7)
colorbar