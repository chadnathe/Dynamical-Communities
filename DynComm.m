function [Clusters] = DynComm(A,kk)

% 
% [Cluster] = DynComm(A,kk) calculates an approximate dynamical cluster 
% partition with adjacency matrix, A, and kk amount of clusters specified. 
% Choose kk from kk = 2...N-1 where N = length(A). 
% 


ns = 100; %choose how many kmedoid samples. Set to 1 for fastest results
N = length(A);
P = [];
Clusters = {1:N};
count = 0;
Cent = cell(1,1);
for k = 2:kk
count = count + 1;
P = zeros(N,k);
for i = 1:N
    for j = 1:length(Clusters)
P(i,j) = sum(A(i,Clusters{j}));
    end
% if k == 2             % If dealing with a graphs which contain isolated
%     if P(i,j) > 0     % nodes, this section may be uncommented.
% P(i,j) = 1;
%     end
% end
end

D = squareform(pdist(P));

    LIST = cell(1,1);
    Stat=zeros(1,ns);
    for ij = 1:ns
        [idx,C,sumd] = kmedoids((1:N)', k, 'Distance', @(ZI, ZJ) D(ZJ, ZI));
        LIST{1,ij} = idx;
        Stat(ij) = sum(sumd);
    end

temp = find(Stat == min(Stat));
CLUS = LIST{1,temp(1)};

K = cell(1,1);
for i = 1:max(CLUS)
a2 = double(CLUS == i);
K{i} = find(a2);
end
Clusters = K;
end


end