% computes k-nearest neighbors for kernels
% data is the are the data values
% query is the (single) point whose neighbors are to be computed
% kernelf is the kernel function
% k is the number of nearest neighbors that have to be computed
% n is the total number of points
function points = kknn(data, query, kernelf, k, n)
%     if(n == k)
%         points = data;
%         return;
%     end
    if(k > n)
        points = horzcat(data, data(:,1:(k - n)));
        return;
    end
    
    dist = zeros(n,1);
    for i = 1:n
        dist(i) = kernelf(data(:,i), data(:,i)) + kernelf(query, query) ...
            - 2*kernelf(data(:,i), query);
    end
    
    [m,ind] = sort(dist);
    points = data(:,ind(1:k));
    
end