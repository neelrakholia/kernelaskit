function covtype_back
clear globals;clc; clear all;
addpath('src/')
addpath('..')

% read train
filename = 'covtype.libsvm.trn.X.bin';
n = 499999;
dim = 54;
train = binread_array(filename, n*dim);
train = reshape(train, dim, n);

% read test
filename = 'covtype.libsvm.tst.X.bin';
m = 80012;
test = binread_array(filename, m*dim);
test = reshape(test, dim, m);

% sample data
n = 2^14;
m = 400;
train = datasample(train, n, 2, 'Replace', false);
test = datasample(test, m, 2, 'Replace', false);

% number of nearest neighbors, and number of trees to generate
K = 10;

% example kernel
sigma = 0.22;

% tree options
maxPointsPerNode = 2^9;
maxLevel        =  12;

% construct tree and search for nearest neighbors
points = zeros(m,K);

tic
root = bsttree_vp(train, 1:n, maxPointsPerNode, maxLevel, sigma, 0, 0);
toc

% search for neighbors for each query point
tic
disteval = 0;
for i = 1:m
    % perform a priority queue based search as described in FLANN
    [point,disev] = psearch(root, train, test(:,i), sigma, K);
    disteval = disteval + disev;
    
    points(i,:) = point; 
end
toc

% actual nearest neighbors using linear search
tic
actual_nn = kknn(train,1:n,test,sigma,K,n);
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
    dist_actual = distk(train(:,actual_nn(i,:)),test(:,i),sigma);
    dist_app = distk(train(:,points(i,:)),test(:,i),sigma);
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