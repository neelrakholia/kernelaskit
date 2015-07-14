function covtype_rand
clear globals;clc; clear all;
addpath('src/')
addpath('..')

% read train
filename = 'covtype.libsvm.trn.X.bin';
n = 2^16;
dim = 54;
train = binread_array(filename, n*dim);
train = reshape(train, dim, n);

% read test
filename = 'covtype.libsvm.tst.X.bin';
m = 6000;
test = binread_array(filename, m*dim);
test = reshape(test, dim, m);

% sample data
n = 2^16;
m = 4000;
train = datasample(train, n, 2, 'Replace', false);
test = datasample(test, m, 2, 'Replace', false);

% number of nearest neighbors, and number of trees to generate
K = 10;
ntree = 20;

% example kernel
sigma = 0.22;

% tree options
maxPointsPerNode = 2^9;
maxLevel        =  12;

% brute force search
tic
actual_nn = kknn(train,1:n,test,sigma,K,n);
toc

% construct tree and search for nearest neighbors
points = zeros(m,K);

% storing previous iteration
test_nn = ones(m, K);

% search for neighbors for each query point
tic
disteval = 0;
treeeval = 0;
k = 0;
acc = 0;
while(k <= ntree && acc < 0.9)
    % construct tree
    root = bsttree_vp(train, 1:n, maxPointsPerNode, maxLevel, sigma, 0, 0);
    for i = 1:m
        % perform tree search
        p = travtree2n(root, test(:,i), sigma);
        
        % search for neighbors
        search_inds = unique([p, test_nn(i,:)]);
        new_nn = kknn(train, search_inds, test(:,i), sigma, K, ...
            numel(search_inds));
        
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
        dist_actual = distk(train(:,actual_nn(i,:)),test(:,i),sigma);
        dist_app = distk(train(:,points(i,:)),test(:,i),sigma);
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