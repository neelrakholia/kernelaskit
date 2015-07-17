classdef bsttree_vp < handle
    
    properties
        data % array of points for our purpose
        ind % indices of data points
        rad = 0 % radius of inner cluster
        cent % point chosen as vantage point
        nsize = 0 % number of elements in the node
        ndepth = 0 % depth of the node
        dise = 0 % number of distance evaluations
        
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
        % diseval: number of distance evaluations
        function root = bsttree_vp(data, indi, msize, mdepth, sigma, depth, diseval)
            % initialize the root
            root.data = data;
            sizep = size(data, 2);
            root.nsize = sizep;
            root.ndepth = depth;
            root.ind = indi;
            
            % check that nsize and ndepth are in acceptable range
            if(root.nsize <= msize || root.ndepth >= mdepth)
                return;
                
                % run kernel k-means and split the data
            else
                % get classification and radius
                [datal, datar, indl, indr, radi, cen, diseval] = classify_vp(data, ...
                    indi, root.nsize, sigma, diseval);
                
                root.dise = diseval;
                
                % assign vantage point
                root.cent = cen;
                
                % assign points
                lchild = datal;
                rchild = datar;
                
                % store radius
                root.rad = radi;
                
                % recursively call bstrree on left and right node
                root.left = bsttree_vp(lchild, indl, msize, mdepth, ...
                    sigma, depth + 1, diseval);
                root.right = bsttree_vp(rchild, indr, msize, mdepth, ...
                    sigma, depth + 1, diseval);
            end % end if
        end % end function
        
    end % end methods
    
    methods
        
        % Priority queue based search for NN. Described in FLANN.
        % root:    the tree containing data
        % data:    all data points
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        % k:       number of nearest neighbors to be found
        function [points,deval] = psearch(root, data, query, sigma, k)
            % initialize dk as infinity
            dk = inf;
            deval = 0;
            
            % call another function to traverse the tree
            [r,dk,deval] = travtree(root, data, query, sigma, k, dk, [], deval);
            
            % select first max number of points
            points = r;
        end
        
        % makes a priority queue based structure while traversing the tree
        % root:    the tree containing data
        % data:    all data points
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        % k:       number of nearest neighbor
        % dk:      distance to furthest neighbor
        % deval:   the number of distance evaluations
        function [q, dk, deval] = travtree(root, data, query, sigma, k, dk, q_in, deval)
            % base case
            % if the root is a leaf
            if(isempty(root.left))
                indi = root.ind;
                q = kknn(data, [indi, q_in], query, sigma, k, ...
                    numel(indi)+numel(q_in));
                deval = deval + numel(indi)+numel(q_in);
                dk_curr = distk(query, data(:,q(k)), sigma);
                if(dk_curr < dk)
                    dk = dk_curr;
                end
                % fprintf('Base case candidate %d, dist: %g\n', q(k), dk_curr);
                return;
                
            end
            
            % get the radius and center
            radius = root.rad;
            center = root.cent;
            
            % calculate distance between query point and center
            dist = distk(query, center, sigma);

            if(dist < radius)
                % store data according to distance from query point
                % fprintf('recursing left: dist %g, radius %g, dk: %g, level %d\n', ...
                    %dist, radius, dk, root.ndepth);
                
                [q, dk, deval] = travtree(root.left, data, query, sigma,...
                    k, dk, q_in, deval);
                
                if(dist + dk > radius)
                    [q, dk, deval] = travtree(root.right, data, query, sigma,...
                        k, dk,q, deval);
                else
                    % fprintf('pruning right: dist %g, radius %g, dk: %g, level %d\n', dist, radius, dk, root.ndepth);
                end
                
            else
                
                % store data according to distance from query point
                % fprintf('recursing right: dist %g, radius %g, dk: %g, level %d\n', dist, radius, dk, root.ndepth);

                [q, dk, deval] = travtree(root.right, data, query, sigma,...
                    k, dk, q_in, deval);
                
                if(dist < radius + dk)
                    [q, dk, deval] = travtree(root.left, data, query, sigma,...
                        k, dk, q, deval);
                else
                    %fprintf('pruning left: dist %g, radius %g, dk: %g, level %d\n', dist, radius, dk, root.ndepth);
                end  
                
                
            end
        end % end function
        
        % travereses to the node closest to the query point
        % root:    the tree containing data
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        function [nn,dev] = travtree2n(root, query, sigma, global_id, data,...
                k, nn, prev, dev)
            % base case
            % if the root is a leaf
            if(isempty(root.left))
                % search for nn
                prev_id = reshape(prev(global_id,:),1,k*numel(global_id));
                search_id = unique([root.ind, prev_id]);
                q = kknn(data, search_id, query, sigma, k, numel(search_id));
                
                % store nn
                nn(global_id,:) = q;
                
                % update computations
                dev = numel(search_id)*numel(global_id);
                return;
            end
            
            % get the radius and center
            radius = root.rad;
            center = root.cent;
            
            % calculate distance between query point and center
            dist = distk(query, center, sigma);
            larr = dist < radius;
            
            % get the right and left queries
            indl = find(larr);
            indr = find(~larr);
            
            % recursive call to whichever center is closer
            if(numel(indl) > 0)
                % store data according to distance from query point
                [nn,dev] = travtree2n(root.left, query(:,indl), sigma, ...
                    global_id(indl), data, k, nn, prev);
            end
            if(numel(indr) > 0)
                % store data according to distance from query point
                [nn,dev] = travtree2n(root.right, query(:, indr), sigma, ...
                    global_id(indr), data, k, nn, prev);
            end % end if
        end % end function
        
        % travereses to the node closest to the query point
        % root:    the tree containing data
        % query:   query point whose NN are to be searched
        % sigma:   bandwidth for rbf kernel
        function q = travtree2n_h(root, query, sigma)
           query = num2cell(query, 1);
           q = cellfun(@(x) travtree2n(root, x, sigma), query, 'UniformOutput', 0);
        end
        
    end % end methods
    
end % end class