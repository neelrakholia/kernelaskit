# kernelaskit
Kernel based distance metrics for ASKIT

I implemented two tree based search algorithms for nearest neighbors in arbitrary feature spaces (kernel distances).
1) Priority Search K-Means Tree Algorithm (described in "Scalable Nearest Neighbor Algorithms for High Dimensional Data" on the FLANN website)
2) Random Sampling (mentioned in "Parallel Algorithms for Nearest Neighbor Search in High Dimensions", a paper written by B. Xiao and G. Biros)

Both algorithms require slightly different tree constructions and hence slight differences in code. Presented below is a brief discussion on each algorithm and the files to execute to test them. 

Priority Search K-Means Tree Algorithm (I used a slight variation of the algorithm because I struggled to implement priority queues in MATLAB)
To run this algorithm we first need a tree, with each node containing data points belonging to a cluster as well as the centers for those clusters. Traditional kernel K-means would not work here because there is no way to determine the geometric center of a cluster in feature space (for instance, the RBF kernel has infinite dimesional feature space). To work around this problem, I used kernel K-medoids to get an approximation of the cluster center. Kernel K-medoids was then called recursively to get a binary search tree. 

To search the tree for nearest neighbors of a query point, the node with the cluster center closest to the point is traversed first. This is recursively performed until a leaf is reached. The points in the leaf node are added to an array. The points in the node adjacent to the leaf node are then added to the array, and so on recrisively. The following code is the gist of the algorithm:

            % base case
            % if the root is a leaf            
            if(isempty(root.left))
               q = root.data;
               return;
            end
            
            % get the cluster centers
            med1 = root.sep(:,1);
            med2 = root.sep(:,2);
            
            % calculate distance between query point and cluster centers
            distl = distk(query, med1, kernelf); % distance to left
            distr = distk(query, med2, kernelf); % distance to right
            
            % recursive call to whichever center is closer 
            if(distl < distr)
                % store data according to distance from query point
                q = horzcat(travtree(root.left, query, kernelf),...
                    root.right.data);
            else
                % store data according to distance from query point
                q = horzcat(travtree(root.right, query, kernelf),...
                    root.left.data);
            end % end if

A fraction of the points in the array are chosen from the top (note that the points are stored in a roughly ascending order of their distance from the query point). Linear K-NN algorithm is run on them and nearest neighbors are calculated. In the extreme case that all points in the array are chosen, we would recover the same nearest neighbors as linear search. 

Random Sampling
The tree construction for this algorithm was a little different. Because we sample for points, there is no need to know the exact cluster centers. As long as we know which points belong to which node, we can construct a binary search tree. Therefore, we can use conventional kernel K-means here to get cluster assignments and get a tree.

To search the tree for nearest neighbors of a query point, we sample ns points from each child node, then visit the node which has the closest sampled point to the query point. This is done recursively until we reach a leaf. The entire process is repeated a number of times to get a set of unique points. Linear knn search is then performed on this set to get the nearest neighbor.

Files for Priority Search K-Means Tree
/src/bsttree_pq.m
/src/distk.m
/src/generate_points.m
/src/kkmeans_m.m
/src/kkmedoids.m
/src/kknn.m

Executable
/test/main1.m

Files for Random Sampling
/src/bsttree_s.m
/src/distk.m
/src/generate_points.m
/src/kkmeans_ap.m
/src/approx_kkmeans.m
/src/kknn.m

Executable
/test/main2.m

Make sure the directory is set to kernelaskit/ before executing the files.

Results
Note that parameters for run can be changed by altering main1.m and main2.m. For the current set of parameters, I get the following values:

Priority Search K-Means Tree
Constructing tree:                         296.698778
Search using Priority Search K-Means Tree: 33.440650
Linear search:                             55.547078  
Accuracy:                                  0.8502

Random Sampling
Constructing tree:                         67.774098
Search using Random Sampling:              42.605046
Linear search:                             53.401546
Accuracy:                                  0.8234

--Accuracy for Priority Search K-Means can be altered by changing the max parameter in main1.m. 
--Accuracy for Random Sampling can be altered by changing the iter and samp parameter in main2.m. 
