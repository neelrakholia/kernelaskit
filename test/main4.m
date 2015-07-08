% Tests for a priority queue based search using vp trees for kernel
% distance metrics
function main4
    clear globals;clc; clear all;  
    addpath('src/')

    % generate points, n = database points, m = query points
    n = 2^14;
    m = 400;
    dim = 8;
    
    % number of nearest neighbors, iterations, and samples to be taken
    K = 10;
    max = floor(0.3*n);
    
    % random generation of database and query points
    point_distribution = 'gaussian';
    r = generate_points(dim, n, point_distribution);
    q = generate_points(dim, m, point_distribution);

    % example kernel
    sigma = 0.5;
    
    % tree options
    maxPointsPerNode = 2^7;
    maxLevel        =  12;
    
    % construct tree and search for nearest neighbors
    points = zeros(m,K);
    
    tic
    root = bsttree_vp(r, 1:n, maxPointsPerNode, maxLevel, sigma, 0);
    toc
    
    % search for neighbors for each query point  
    tic
    for i = 1:m
        % perform a priority queue based search as described in FLANN
        point = psearch(root, q(:,i), sigma, max);
        
        % search for nearest neighbors among those points
        len = size(point);
        len = len(2);
        points(i,:) = kknn(r, point, q(:,i), sigma, K, len);
    end
    toc
    
    % actual nearest neighbors using linear search
    tic
    actual_nn = kknn(r,1:n,q,sigma,K,n);
    toc
    
    distk(q(:,1),r(:,actual_nn(1,:)), sigma)
    
    suml = 0;
    % calculate accuracy by comparing neighbors found
    for i = 1:m
        suml = suml + length(intersect(points(i,:), actual_nn(i,:)));
    end
    
    % print accuracy
    fprintf('Accuracy: %f\n', suml/(m*K))
    
    % display fraction of distance evaluations conducted
    fprintf('Distance evaluations: %f\n', (max*m)/(m*n))

end