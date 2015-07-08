% computes k-nearest neighbors for kernels
% data:    data values
% query:   point whose neighbors are to be computed
% sigma:   bandwidth
% k:       number of nearest neighbors that have to be computed
% n:       total number of points
function points = kknn(data, indi, query, sigma, k, n)
% in the event that number of neighbors is greater than n
if(k > n)
    points = horzcat(data, data(:,1:(k - n)));
    return;
end

% linear search through all the points
dist = distk(data(:, indi), query, sigma);

len = size(query, 2);

% if there is 1 query point or many
if(len == 1)
    % sort distance and report k NN
    [m,ind] = sort(dist);
    points = indi(ind(1:k));
else
    % sort distance and report k NN
    [m,ind] = sort(dist);
    ind = ind';
    % datasum = sum(data, 1);
    points = indi(ind(:,1:k));
end

end
