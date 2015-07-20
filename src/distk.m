% computes the distance for kernels
% p1:      point1
% p2:      point2
% kernelf: kernel function used in calculating distanced
function d = distk(X, x, sigma)
% select kernel type
kerneltype = 'rbf';

% calculate distance
if(strcmp(kerneltype,'rbf'))
    d = rbf(X, x, sigma);
else
    d = polyk(X, x, sigma);
end

end