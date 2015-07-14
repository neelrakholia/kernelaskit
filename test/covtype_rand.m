function covtype_rand
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
m = 81012;
test = binread_array(filename, m*dim);
test = reshape(test, dim, m);


% number of nearest neighbors, and number of trees to generate
K = 10;
ntree = 20;

% example kernel
sigma = 0.3;

tic
actual_nn = kknn(train,1:n,test,sigma,K,n);
toc

% tree options
maxPointsPerNode = 2^13;
maxLevel        =  12;

% construct tree and search for nearest neighbors
points = zeros(m,K);
point = zeros(m, K*ntree);

% actual nearest neighbors using linear search
tic
actual_nn = kknn(train,1:n,test,sigma,K,n);
toc

% storing previous iteration
test_nn = ones(m, K);

% search for neighbors for each query point
tic
disteval = 0;
treeeval = 0;
k = 0;
acc = 0;
while(k <= ntree && acc ~= 1.0)
    % construct tree
    root = bsttree_vp(r, 1:n, maxPointsPerNode, maxLevel, sigma, 0, 0);
    for i = 1:m
        % perform tree search
        p = travtree2n(root, q(:,i), sigma);
        
        % search for neighbors
        search_inds = unique([p, test_nn(i,:)]);
        new_nn = kknn(r, search_inds, q(:,i), sigma, K, numel(search_inds));
        
        % update disteval
        disteval = disteval + numel(search_inds);
        
        % update array to store current iteration
        test_nn(i,:) = new_nn;
        
        % store points
        points(i,:) = test_nn(i,:);
    end
    
    treeeval = treeeval + root.dise;
    
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
    acc = suml/(m*K);
    fprintf('Accuracy: %f\n', suml/(m*K));
    
    % display fraction of distance evaluations conducted while constructing
    % the tree
    fprintf('Distance evaluations in tree: %f\n', treeeval/(m*n));
    
    % display fraction of distance evaluations conducted
    fprintf('Distance evaluations: %f\n', disteval/(m*n));
    
    % display total number of distance evaluations:
    fprintf('Total distance evaluations: %f\n', disteval/(m*n) + ...
        treeeval/(m*n));
    
    % print ratio of distance
    fprintf('Average ratio of distance: %f\n',mean(distr));
    
    k = k + 1;
end
toc