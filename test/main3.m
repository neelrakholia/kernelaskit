% Tests for a random sampling based nearest neighbor search for kernel
% distance metrics. Here we use Naive Trees. 
function main3
    clear globals;clc; clear all;  
    addpath('src/')

    % generate points, n = database points, m = query points
    n = 2^12;
    m = 400;
    dim = 8;
    
    % number of nearest neighbors, iterations, and samples to be taken
    K = 10;
    iter = 25;
    samp = 2^3;
    
    % random generation of database and query points
    point_distribution = 'gaussian';
    r = generate_points(dim, n, point_distribution);
    q = generate_points(dim, m, point_distribution);

    % example kernel
    sigma = 2;
    
    % tree options
    maxPointsPerNode = 2^7;
    maxLevel        =  10;
    
    % construct tree and search for nearest neighbors
    points = zeros(m,K);
    
    tic
    root = bsttree_s(r, maxPointsPerNode, maxLevel, sigma, 0, 1);
    toc
    
    % search for neighbors for each query point
    disteval = 0;
    tic
    for i = 1:m
        % iterate through the tree by sampling different points everytime
        point = root.sampsearch(q(:,i), sigma, K, samp);
        for j = 1:iter - 1
            point = horzcat(point, root.sampsearch(q(:,i), sigma, K, samp));
        end
        
        % search for nearest neighbors among gathered points
        point = unique(point.', 'rows');
        point = point.';
        len = size(point);
        len = len(2);
        disteval = disteval + len;
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
    fprintf('Distance evaluations: %f\n', disteval/(m*n))

end