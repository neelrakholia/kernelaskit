% Runs kernel k-medoid and outputs labels and centroids
% points:  data to be clustered
% kernelf: kernel function
% k:       number of clusters to be formed
% n:       total number of points
function [label, medoid] = kkmeans_m(points, kernelf, k, n)

% get distance matrix
K = zeros(n, n);
for i = 1:n
    for j = 1:n
        K(i, j) = distk(points(:, i), points(:, j), kernelf);
    end
end

% run kernel k-medoid
cluster = kkmedoid(points', k, K);
label = cluster.label;
medoid = cluster.medoids;

end