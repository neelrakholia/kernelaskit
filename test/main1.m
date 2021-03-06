% Tests for a priority queue based nearest neighbor search for kernel
% distance metrics
function main1
    clear globals; clc; clear all;  
    addpath('src/')

    % generate points, n = database points, m = query points
    n = 2^12;
    m = 400;
    dim = 8;
    
    % number of nearest neighbors and maximum number of points to scan
    K = 10;
    max = floor(0.75*n);
    
    % random generation of database and query points
    point_distribution = 'gaussian';
    r = generate_points(dim, n, point_distribution);
    q = generate_points(dim, m, point_distribution);

    % example kernel
    sigma = 2;
    kernelff = @(x, y) exp(-norm(x - y, 2)^2/(2*sigma^2));
    
    % tree options
    maxPointsPerNode = 2^7;
    maxLevel        =  10;
    
    % construct tree and search for nearest neighbors
    points = zeros(m,K);
    
    tic
    root = bsttree_pq(r, maxPointsPerNode, maxLevel, kernelff, 0);
    toc
    
    % search for neighbors for each query point
    tic
    for i = 1:m
        % perform a priority queue based search as described in FLANN
        point = psearch(root, q(:,i), sigma, max);
        
        % search for nearest neighbors among those points
        len = size(point);
        len = len(2);
        points(i,:) = sum(kknn(point, q(:,i), sigma, K, len))';
    end
    toc
    
    % actual nearest neighbors using linear search
    tic
    actual_nn = kknn(r,q,sigma,K,n);
    toc
    
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