% Tests for a node based search using vp trees for kernel
% distance metrics
function main5
    clear globals;clc; clear all;  
    addpath('src/')

    % generate points, n = database points, m = query points
    n = 2^14;
    m = 400;
    dim = 8;
    
    % number of nearest neighbors, and number of trees to generate
    K = 10;
    ntree = 150;
    
    % random generation of database and query points
    point_distribution = 'gaussian';
    r = generate_points(dim, n, point_distribution);
    q = generate_points(dim, m, point_distribution);

    % example kernel
    sigma = 2;
    
    % tree options
    maxPointsPerNode = 2^7;
    maxLevel        =  12;
    
    % construct tree and search for nearest neighbors
    points = zeros(m,K);
    
    % construct a number of trees
    root(1) = bsttree_vp(r, maxPointsPerNode, maxLevel, sigma, 0); 
    tic
    for j = 2:ntree
    root(j) = bsttree_vp(r, maxPointsPerNode, maxLevel, sigma, 0);
    end
    toc
    
    % search for neighbors for each query point  
    tic
    disteval = 0;
    for i = 1:m
        % perform tree search
        point = travtree2n(root(1), q(:,i), sigma);
        for k = 2:ntree
            point = horzcat(point, travtree2n(root(k), q(:,i), sigma));
        end
        
        % search for nearest neighbors among gathered points
        point = unique(point.', 'rows');
        point = point.';
        
        % search for nearest neighbors among those points
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