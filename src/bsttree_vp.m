classdef bsttree_vp < handle
    
    properties
        data % array of points for our purpose
        ind % indices of data points
        rad = 0 % radius of inner cluster
        cent % point chosen as vantage point
        nsize = 0 % number of elements in the node
        ndepth = 0 % depth of the node
        
    end
    
    properties  (SetAccess = private)
        left = bsttree_vp.empty; % left subtree
        right = bsttree_vp.empty; % right subtree
    end
    
    methods
        
        % initialization
        % data:    data from which the tree needs to be constructed
        % msize:   maximum number of data points allowed in a node
        % mdepth:  maximum depth of the tree
        % kernelf: rbf kernel bandwidth
        % depth:   intial depth of the tree
        function root = bsttree_vp(data, indi, msize, mdepth, sigma, depth)
            % initialize the root
            root.data = data;
            sizep = size(data);
            sizep = sizep(2);
            root.nsize = sizep;
            root.ndepth = depth;
            root.ind = indi;
            
            % check that nsize and ndepth are in acceptable range
            if(root.nsize <= msize || root.ndepth >= mdepth)
                return;
                
                % run kernel k-means and split the data
            else
                % get classification and radius
                [datal, datar, indl, indr, radi, cen] = classify_vp(data, ...
                    indi, root.nsize, sigma); 
                
                % assign vantage point
                root.cent = cen;
                
                % assign points
                lchild = datal;
                rchild = datar;
                
                % store radius
                root.rad = radi;
                
                % recursively call bstrree on left and right node
                root.left = bsttree_vp(lchild, indl, msize, mdepth, ...
                    sigma, depth + 1);
                root.right = bsttree_vp(rchild, indr, msize, mdepth, ...
                    sigma, depth + 1);
            end % end if
        end % end function
        
    end % end methods
    
    methods
        
        % Priority queue based search for NN. Described in FLANN.
        % root:    the tree containing data
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        function points = psearch(root, query, sigma, max)           
           % call another function to traverse the tree
           r = travtree(root, query, sigma);
           
           % select first max number of points
           points = r(:,1:max);        
        end
        
        % makes a priority queue based structure while traversing the tree
        % root:    the tree containing data
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        function q = travtree(root, query, sigma)
            % base case
            % if the root is a leaf            
            if(isempty(root.left))
               q = root.ind;
               return;
            end
            
            % get the radius and center
            radius = root.rad;
            center = root.cent;
            
            % calculate distance between query point and center
            dist = distk(query, center, sigma);
            
            % recursive call to whichever center is closer 
            if(dist < radius)
                % store data according to distance from query point
                q = horzcat(travtree(root.left, query, sigma),...
                    root.right.ind);
            else
                % store data according to distance from query point
                q = horzcat(travtree(root.right, query, sigma),...
                    root.left.ind);
            end % end if
        end % end function
        
        % travereses to the node closest to the query point
        % root:    the tree containing data
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        function q = travtree2n(root, query, sigma)
            % base case
            % if the root is a leaf            
            if(isempty(root.left))
               q = root.ind;
               return;
            end
            
            % get the radius and center
            radius = root.rad;
            center = root.cent;
            
            % calculate distance between query point and center
            dist = distk(query, center, sigma);
            
            % recursive call to whichever center is closer 
            if(dist < radius)
                % store data according to distance from query point
                q = travtree2n(root.left, query, sigma);
            else
                % store data according to distance from query point
                q = travtree2n(root.right, query, sigma);
            end % end if
        end % end function

    end % end methods
  
end % end class