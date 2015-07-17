function covtype_rand
clear globals; clc; clear all;
addpath('src/')
addpath('..')

% read train
filename = 'susy.icml.trn.X.bin';
n = 4499999;
dim = 18;
train = binread_array(filename, n*dim);
train = reshape(train, dim, n);

% read test
filename = 'susy.icml.tst.X.bin';
m = 499999;
test = binread_array(filename, m*dim);
test = reshape(test, dim, m);

% sample data
n = 4400000;
m = 40;
train = datasample(train, n, 2, 'Replace', false);
test = datasample(test, m, 2, 'Replace', false);

% number of nearest neighbors, and number of trees to generate
K = 10;
ntree = 20;

% example kernel
sigma = 0.15;

% tree options
maxPointsPerNode = 2^13;
maxLevel        =  12;

% brute force search
tic
piece = 10;
actual_nn = kknn(train,1:n,test(:,1:m/piece),sigma,K,n);
for i = 2:piece
    actual_nn = vertcat(actual_nn, kknn(train,1:n,...
        test(:,(i - 1)*(m/piece) + 1:i*m/piece),sigma,K,n));
end
toc

% construct tree and search for nearest neighbors
points = zeros(m,K);

% storing previous iteration
test_nn = ones(m, K);

% array for storing time
elapsed_time_array = [];

% search for neighbors for each query point
disteval = 0;
treeeval = 0;
k = 0;
acc = 0;

% iterate through all the the trees
while(k <= ntree && acc < 0.9)
    % construct tee
    tic
    root = bsttree_vp(train, 1:n, maxPointsPerNode, maxLevel, sigma, 0, 0);
    elapsed_time_array(end + 1) = toc;
    
    % search tree
    tic
    [new_nn,deval] = travtree2n(root, test, sigma, ...
        1:m, train, K, points, test_nn, 0);
    test_nn = new_nn;
    points = test_nn;
    disteval = disteval + deval;
    elapsed_time_array(end + 1) = toc;
    
    % evaluate performace
    tic
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
    elapsed_time_array(end + 1) = toc;
end

fprintf('Time spent constructing trees: %g\n',...
    sum(elapsed_time_array(1:3:end)));
fprintf('Time spent searching trees: %g\n',...
    sum(elapsed_time_array(2:3:end)));
fprintf('Total time spent running algo: %g\n',...
    sum(elapsed_time_array(1:3:end)) + sum(elapsed_time_array(2:3:end)));
fprintf('Total time spent running: %g\n',...
    sum(elapsed_time_array))