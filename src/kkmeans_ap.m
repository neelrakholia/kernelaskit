% Runs kernel k-means and outputs labels
% points:  data to be clustered
% kernelf: kernel function
% k:       number of clusters to be formed
% n:       total number of points
function label = kkmeans_ap(points, kernelf, k, n)
% set maimum number of iterations
max_iter = 100;

% the number of sample points that we take
m = floor(n/4);

% sample m data points
perm = randperm(n);
indices = perm(1:m);

% get kernel matrix
K = zeros(m, n);
for i = 1:m
    for j = 1:n
        K(i, j) = kernelf(points(:, indices(i)), points(:, j));
    end
end

% run approximate kernel k-means algorithm
label = approx_kkmeans(K, k, max_iter, indices);
end

