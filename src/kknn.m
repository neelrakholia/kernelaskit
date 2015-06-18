% computes k-nearest neighbors for kernels
% data:    data values
% query:   point whose neighbors are to be computed
% kernelf: kernel function
% k:       number of nearest neighbors that have to be computed
% n:       total number of points
function points = kknn(data, query, kernelf, k, n)
    % in the event that number of neighbors is greater than n
    if(k > n)
        points = horzcat(data, data(:,1:(k - n)));
        return;
    end
    
    % linear search through all the points
    dist = zeros(n,1);
    for i = 1:n
        dist(i) = distk(data(:,i), query, kernelf);
    end
    
    % sort distance and report k NN
    [m,ind] = sort(dist);
    points = data(:,ind(1:k));
    
end