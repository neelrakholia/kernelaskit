% computes the distance for kernels
% p1:      point1
% p2:      point2
% kernelf: kernel function used in calculating distanced
function d = distk(X, x, sigma)

% evaluate euclidean part
[D,N] = size(X);
[D,n] = size(x);
X2 = sum(X.^2,1);
x2 = sum(x.^2,1);
dotProd = X'*x;
d = repmat(x2,N,1)+repmat(X2',1,n)-2*dotProd;

% evaluate kernel
% evaluate kernel
d = exp(-d/(2*sigma^2));
d = sqrt(-2*d + 2);

end