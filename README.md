### kernelaskit
### Author:     Neel Rakholia
### Co-Author:  Bill March
### Institute:   ICES, UT Austin, Austin, Texas

I implemented two tree based search algorithms for nearest neighbors in 
arbitrary feature spaces (kernel distances).

1) Priority Search K-Means Tree Algorithm (described in "Scalable Nearest 
Neighbor Algorithms for High Dimensional Data" on the FLANN website)
2) Random Sampling (mentioned in "Parallel Algorithms for Nearest Neighbor 
Search in High Dimensions", a paper written by B. Xiao and G. Biros)

Both algorithms require slightly different tree constructions and hence slight 
differences in code. Presented below is a brief discussion on each algorithm 
and the files to execute to test them. 


Priority Search K-Means Tree Algorithm (I used a slight variation of the 
algorithm because I struggled to implement priority queues in MATLAB)

To run this algorithm we first need a tree, with each node containing data 
points belonging to a cluster as well as the centers for those clusters. 
Traditional kernel K-means would not work here because there is no way to 
determine the geometric center of a cluster in feature space (for instance, 
the RBF kernel has infinite dimesional feature space). To work around this 
problem, I used kernel K-medoids to get an approximation of the cluster center.
 Kernel K-medoids was then called recursively to get a binary search tree. 

To search the tree for nearest neighbors of a query point, the node with the 
cluster center closest to the point is traversed first. This is recursively 
performed until a leaf is reached. The points in the leaf node are added to
 an array. The points in the node adjacent to the leaf node are then added to
 the array, and so on recrisively. The following code is the gist of the algorithm:

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

A fraction of the points in the array are chosen from the top (note that the 
points are stored in a roughly ascending order of their distance from the query
point). Linear K-NN algorithm is run on them and nearest neighbors are calculated.
In the extreme case that all points in the array are chosen, we would recover the 
same nearest neighbors as linear search. 


Random Sampling

The tree construction for this algorithm was a little different. Because we 
sample for points, there is no need to know the exact cluster centers. As long 
as we know which points belong to which node, we can construct a binary search 
tree. Therefore, we can use conventional kernel K-means here to get cluster 
assignments and get a tree.

To search the tree for nearest neighbors of a query point, we sample ns points 
from each child node, then visit the node which has the closest sampled point 
to the query point. This is done recursively until we reach a leaf. The entire 
process is repeated a number of times to get a set of unique points. Linear knn
 search is then performed on this set to get the nearest neighbor.


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

Note that parameters for run can be changed by altering main1.m and main2.m. 
For the current set of parameters, I get the following values:

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

--Accuracy for Priority Search K-Means can be altered by changing the max 
parameter in main1.m. 
--Accuracy for Random Sampling can be altered by changing the iter and samp 
parameter in main2.m. 

Note: To achieve > 0.9 accuracy, we need to use 80 percent of the points for 
Priority Search, and 35 iterations for Random Sampling


UPDATE 1

I implemented a 3rd algorithm. It involves "naive" construction of trees. What
this essentially means is that I select two random points from a node and make 
them cluster centers. The points closest to a cluster center are assigned to 
that cluster. This obviously results in a highly unbalanced tree structure, but
the construction is extremely fast. The search algorithm for this tree structure
is the same as Random Sampling. 

Executable: /test/main3.m


Results

Constructing tree:                         1.805226
Search using Random Sampling:              43.837278
Linear search:                             53.887611 
Accuracy:                                  0.6455


UPDATE 2

Optimized linear search for algorithms by removing for loops. Increased the 
speed by more than 100x. Also implemented an additional algorithm for tree 
construction (vp-trees). 


Vantage Point Trees (vp-trees)

A vantage-point tree, or VP tree is a BSP tree that segregates data in a metric
space by choosing a position in the space (the "vantage point") and dividing the
data points into two partitions: those that are nearer to the vantage point than
a threshold, and those that are not. By repeatedly applying this procedure to 
partition the data into smaller and smaller sets, a tree data structure is 
created where neighbors in the tree are likely to be neighbors in the space. 
(Wikipedia) I used kernel distances for tree construction and a priority queue 
based search for traversing the tree. The algorithm for selecting the vantage 
points is described in the paper: "Data Structures and Algorithms for Nearest 
Neighbor Search in General Metric Spaces" by P. Yianilos. Essentially the data 
point that maximizes the variance in distance is selected.


Files for Vantage Point Trees
/src/bsttree_vp.m
/src/distk.m
/src/generate_points.m
/src/classify_vp.m
/src/kknn.m

Executable
/test/main4.m

The updated numbers for all algorithms are presented below:


Results

Priority Search K-Means Tree
Constructing tree:                         295.370181
Search using Priority Search K-Means Tree: 0.531603
Linear search:                             0.080025 
Accuracy:                                  0.8818

Random Sampling
Constructing tree:                         72.079934
Search using Random Sampling:              28.777789
Linear search:                             0.080642
Accuracy:                                  0.8153

"Naive" Tree Construction and Search
Constructing tree:                         1.875423
Search using Random Sampling:              25.557008
Linear search:                             0.070205
Accuracy:                                  0.6253

VP-Trees 
Constructing tree:                         0.520415
Search using Random Sampling:              1.048315
Linear search:                             0.324681
Accuracy:                                  0.9728
(note these times are for 2^14 points, 4 times the number of points used for testing other algorithms)

VP-Trees seem promising at this point because tree construction time, search 
time, and even accuracy is the best for them among all tested algorithms. 


UPDATE 3

I discovered a couple of bugs in the file distk.m and classify_vp.m. These were
first corrected. The way I calculated accuracy was also ineffecient and faulty.
To overcome this, I implemented a global indexing scheme for the data points. 
With VP-trees giving the highest accuracy and fastest construction time, 
I decided to explore them further. In particular, I decided to try two seperate
algoithms:

1) Backtrack search: Searches through the tree until it finds the exact neighbor.
Guarantees an accuracy of 1.
2) Random Tree search: A number of VP-trees are constructed and for each query 
a greedy search is performed on each tree. The larger the number of trees used,
the greater is the number of actual neighbors discovered. 

To test how well Random Tree search performs, I also calculated the average of 
the ratio of distance between app. NN and actual NN. In addition to this, I also
calculated what fraction of total distance calculations (performed by linear 
brute force search algorithm) are computed by both the aforementioned algorithms.
The results are presented below:


Files for Backtrack Search
/src/bsttree_vp.m
/src/distk.m
/src/generate_points.m
/src/classify_vp.m
/src/kknn.m

Executable
/test/main4.m

Files for Random Tree Search
/src/bsttree_vp.m
/src/distk.m
/src/generate_points.m
/src/classify_vp.m
/src/kknn.m

Executable
/test/main5.m


Results

Backtrack Search
Constructing tree:                         0.167403
Search using backtrack search:             16.510435
Linear search:                             0.362102 
Accuracy:                                  1.0000
Distance evaluations in tree:              0.002500
Distance evaluations:                      1.064754
Total distance evaluations:                1.067254

Random tree search (20 trees)
Constructing trees + search:               13.046991
Linear search:                             0.362944
Accuracy:                                  0.9403
Distance evaluations in tree:              0.050000
Distance evaluations:                      0.165302
Total distance evaluations:                0.215302
Average ratio of distance:                 1.004181

As expected there is a trade-off between accuracy and number of distance 
evaluations. Backtrack search for high dimensions is inevitably slow because of
the curse of dimensionality. Random tree search is faster but not completely 
accurate. As we study other data structures, the values listed above will be 
used as a standard for comparison. 


UPDATE 4

I decided to experiment a little with the parameters, and see how that affects 
the performance of both the methods described in UPDATE 3. In particular, I
changed the number of dimensions and the RBF kernel bandwidth, sigma. Everything
else was kept constant. The results are presented below:

http://bit.ly/1UUP5UA

