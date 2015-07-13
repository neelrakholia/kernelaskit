function realdata_back
clear globals;clc; clear all;
addpath('src/')

% generate points, n = database points, m = query points
n = 2^14;
m = 400;
dim = 4;
dim_art = 100;

% number of nearest neighbors, iterations, and samples to be taken
K = 10;

% random generation of database and query points
% random generation of database and query points
point_distribution = 'gaussian';
r = generate_points(dim, n, point_distribution);
q = generate_points(dim, m, point_distribution);
pad = zeros(dim_art - dim, n);
pad_q = zeros(dim_art - dim, m);
[rmatrix,R] = qr(randn(dim_art));
rmatrix = rmatrix*diag(sign(diag(R)));
r = vertcat(r, pad);
r = rmatrix * r;
q = vertcat(q, pad_q);
q = rmatrix * q;

% example kernel
sigma = 2;

% tree options
maxPointsPerNode = 2^7;
maxLevel        =  12;

% construct tree and search for nearest neighbors
points = zeros(m,K);

tic
root = bsttree_vp(r, 1:n, maxPointsPerNode, maxLevel, sigma, 0, 0);
toc

% search for neighbors for each query point
tic
disteval = 0;
for i = 1:m
    % perform a priority queue based search as described in FLANN
    [point,disev] = psearch(root, r, q(:,i), sigma, K);
    disteval = disteval + disev;
    
    % search for nearest neighbors among those points
    %         len = size(point);
    %         len = len(2);
    points(i,:) = point; %kknn(r, point, q(:,i), sigma, K, len);
end
toc

% actual nearest neighbors using linear search
tic
actual_nn = kknn(r,1:n,q,sigma,K,n);
toc

suml = 0;
% calculate accuracy by comparing neighbors found
if(K == 1)
    for i = 1:m
        suml = suml + length(intersect(points(i), actual_nn(i)));
    end
else
    for i = 1:m
        suml = suml + length(intersect(points(i,:), actual_nn(i,:)));
    end
end

% calculate distance ratio
distr = zeros(m, 1);
for i = 1:m
    dist_actual = distk(r(:,actual_nn(i,:)),q(:,i),sigma);
    dist_app = distk(r(:,points(i,:)),q(:,i),sigma);
    distr(i) = mean(dist_app ./ dist_actual);
end

% print accuracy
fprintf('Accuracy: %f\n', suml/(m*K))

% display fraction of distance evaluations conducted while constructing
% tree
fprintf('Distance evaluations in tree: %f\n', root.dise/(m*n))

% display fraction of distance evaluations conducted
fprintf('Distance evaluations: %f\n', disteval/(m*n))

% display total number of distance evaluations:
fprintf('Total distance evaluations: %f\n', disteval/(m*n) + ...
    root.dise/(m*n));

% print ratio of distance
fprintf('Average ratio of distance: %f\n',mean(distr));

end