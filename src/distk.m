% computes the distance for kernels
% p1:      point1
% p2:      point2
% kernelf: kernel function used in calculating distanced
function d = distk(p1, p2, kernelf)
d = kernelf(p1, p1) + kernelf(p2, p2) - 2*kernelf(p1, p2);
d = d/(d + 1);
end