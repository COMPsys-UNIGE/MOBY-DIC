function object = computeTree(object,varargin)
% computeTree   Computes the binary search tree associated to the polytopic partition
%
% OBJ = computeTree(OBJ)
% Computes the bynary search tree with the default settings (see below).
% OBJ is the pwagFunction object.
%
% OBJ = computeTree(OBJ,MODE)
% Computes the binary search tree by defining a MODE which can be either
% 'regions' (default) or 'functions'.
% In the first case the domain is split based on the polytopes, in the
% second case the domain is split based on the affine functions to be
% evaluated. This means that if there are only two regions on a side of an
% edge and in these two regions the affine functions to be computed are the
% same, the two regions are considered as only one region.
%
% OBJ = computeTree(OBJ,BALANCE)
% Computes the binary search tree by specifying a value which determines
% the balance of the tree itself. BALANCE is a number between 0 and 1
% which indicates if you prefer to have a more balanced tree (balance -> 1)
% but with a greater number of nodes or a less balanced tree (balance -> 0)
% with less nodes. Default value 0.5.
%
% OBJ = computeTree(OBJ,MODE,BALANCE)
% Computes the binary search tree by specifying both mode and balance.
%
% See also: pwagFunction/plotTree.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA  02111-1307  USA


if nargin < 1 || nargin > 3
    error('Wrong input arguments');
end

if nargin == 1
    mode = 'regions';
    balance = 0.5;
elseif nargin == 2
    if ischar(varargin{1})
        mode = varargin{1};
        if ~strcmpi(mode,'regions') && ~strcmpi(mode,'functions')
            error('mode must be either ''regions'' or ''functions''');
        end
        balance = 0.5;
    elseif varargin{1} >= 0 && varargin{1} <= 1
        balance = varargin{1};
        mode = 'regions';
    else
        error('Input argument must be a string or a number between 0 and 1');
    end
else
    mode = varargin{1};
    if ~strcmpi(mode,'regions') && ~strcmpi(mode,'functions')
        error('mode must be either ''regions'' or ''functions''');
    end
    balance = varargin{2};
    if balance < 0 || balance > 1
        error('balance must be a number between 0 and 1');
    end
end

% All available edges
HH = object.edges.H;
KK = object.edges.K;

% Domain edges
[Hd, Kd] = object.getDomain();

% Create edge set
edge_set = 1:size(HH,1);

MOBYDICpars = getMOBYDICpars;
tol = MOBYDICpars.roundtol;

% Remove domain edges from edge set
HK = [HH KK];
HKd = [Hd Kd];
HK = tol*round(HK/tol);
HKd = tol*round(HKd/tol);
rem = ismember(HK,HKd,'rows');
edge_set(rem) = [];

% All the regions
regions = object.regions;

% Number of regions
nr = numel(regions);

% Root node initialization
nodes(1).name = 1;
nodes(1).edge_index = [];
nodes(1).region_index = 1:nr;
nodes(1).edge_set = edge_set;
nodes(1).children = [2 3];
nodes(1).leaf = 0;
nodes(1).depth = 0;

% Indexes of all the nodes (it is updated at each iteration)
nodeset = 1;

% Index of the current element in nodeset
inodeset = 1;

% Array containing the depth of each node
depths = [];

% When proceed becomes 0 the tree has been completely computed
proceed = 1;

% Tolerance on the Chebyshev radius. A region with Chebyshev radius lower
% than tol is considered empty
MOBYDICpars = getMOBYDICpars;
tol = MOBYDICpars.emptytol;

% Main loop, constructing all nodes

h = waitbar(0,'Computing binary search tree...');

added_edges = [];
nedges = numel(edge_set);

while proceed
    
    % Creating ii-th node
    ii = nodeset(inodeset);
    
    %     disp([num2str(ii),' - ',num2str(nodeset(end))])
    
    % If the current node is not a leaf...
    if nodes(ii).leaf == 0
        edge_set = nodes(ii).edge_set;
        region_index = nodes(ii).region_index;
        
        cost = zeros(numel(region_index),numel(edge_set));
        cx = cost;
        
        % Find the edge used to split the domain in 2 parts.
        % The edge is chosen to minimize this cost: bal*|sx-dx|+(1-bal)*cx
        % where sx is the number of regions in the "left" of the current edge,
        % dx the number of regions in the right and cx the number of regions crossed
        % by the edge. This formula is chosen in order to have a "balanced" tree.
        
        for i = 1:numel(edge_set)
            
            % Current edge
            Hedge = HH(edge_set(i),:);
            Kedge = KK(edge_set(i));
            
            if strcmpi(mode,'functions')
                % Looping in all regions which can still be considered (the ones present
                % in region_index)
                fsx = zeros(0,object.getCodomainDimensions);
                fdx = zeros(0,object.getCodomainDimensions);
                for j = region_index
                    
                    % Edges of the current region
                    [Hloc, Kloc] = object.getEdges(j);
                    
                    % Chebyshev radius of the part of the current region which
                    % is in the "left" of the current edge
                    rsx = chebyradius([Hloc;Hedge],[Kloc;Kedge]);
                    
                    % If the radius of the partial region is greater than tol,
                    % sx is set to 1 (this means that there is a significant 
                    % portion of the region at the left of the current edge)
                    sx = 0;
                    if rsx > tol
                        if ~ismember(regions(j).Ifunctions',fsx,'rows')
                            fsx = [fsx; regions(j).Ifunctions']; 
                            sx = 1;
                        end
                    end
                    
                    % The same for the "right" side
                    rdx = chebyradius([Hloc;-Hedge],[Kloc;-Kedge]);
                    dx = 0;
                    if rdx > tol
                        if ~ismember(regions(j).Ifunctions',fdx,'rows')
                            fdx = [fdx; regions(j).Ifunctions']; 
                            dx = 1;
                        end
                    end                
                    
                    % Cost(j,i) is 1 if the region is totally at the "left" side of the
                    % current edge, -1 if it is totally in the right side, 0 if the edge
                    % crosses
                    cost(j,i) = sx-dx;
                    cx(j,i) = sx+dx == 2;
                end
                
            else
                
                % Looping in all regions which can still be considered (the ones present
                % in region_index)
                for j = region_index
                    
                    % Edges of the current region
                    [Hloc, Kloc] = object.getEdges(j);
                    
                    % Chebyshev radius of the part of the current region which
                    % is in the "left" of the current edge
                    rsx = chebyradius([Hloc;Hedge],[Kloc;Kedge]);
                    
                    % If the radius of the partial region is greater than a
                    % tolerance, sx is set to 1 (this means that there is
                    % a significant portion of the region at the left of the current edge)
                    sx = rsx > tol;
                    
                    % The same for the "right" side
                    rdx = chebyradius([Hloc;-Hedge],[Kloc;-Kedge]);
                    dx = rdx > tol;
                    
                    % Cost(j,i) is 1 if the region is totally at the "left" side of the
                    % current edge, -1 if it is totally in the right side, 0 if the edge
                    % crosses
                    cost(j,i) = sx-dx;
                    
                    % cx(j,i) is 1 if edge i splits region j into two parts
                    cx(j,i) = sx+dx == 2;
                end
            end
        end
        
        % The cost associated to each edge is defined as bal*|sx-dx|+(1-bal)cx
        tmp = sum(cx);
        cost = sum(cost);
        cost = balance*abs(cost) + (1-balance)*tmp;
        
        % Idx contains the index of the edge with smaller cost
        [min_cost, idx] = min(cost); %#ok<ASGLU>
        
        % The index with smaller cost is chosen to split the domain
        edge_index = nodes(ii).edge_set(idx);
        
        added_edges = unique([added_edges edge_index]);
        waitbar(numel(added_edges)/nedges);
        
        % Chosen edge
        Hedge = HH(edge_index,:);
        Kedge = KK(edge_index);
        
        % Indexes of the current node's children
        child1 = nodes(ii).children(1);
        child2 = nodes(ii).children(2);
        
        region_index1 = [];
        region_index2 = [];
        edge_set1 = [];
        edge_set2 = [];
        
        % Looping in all regions which can still be considered (the ones present
        % in region_index)
        for j = nodes(ii).region_index
            
            % Edges of the current region
            [Hloc, Kloc] = object.getEdges(j);
                       
            % Chebyshev radius of the parts of the current region which
            % are in the "left" and "right" of the chosen edge
            rsx = chebyradius([Hloc;Hedge],[Kloc;Kedge]);
            rdx = chebyradius([Hloc;-Hedge],[Kloc;-Kedge]);
            
            % If the part of the region in the "left" or "right" is too small
            % with respect to a given tolerance, it is marked
            % to be removed
            emptysx = rsx < tol;
            emptydx = rdx < tol;
            
            % Edges to add in edge_set
            edges_to_add = regions(j).Iedges(:,1)';
            
            % Remove edges that are not edges of the parent
            rem = ~ismember(edges_to_add,nodes(ii).edge_set);
            edges_to_add(rem) = [];
            
            % If the part of the region in the "left" of the chosen edge is not
            % too small, add the edges of the current region to the edge_set field of
            % the second child of the current node
            if ~emptysx
                region_index2 = [region_index2 j]; %#ok<*AGROW>
                edge_set2 = [edge_set2 edges_to_add];
            end
            % If the part of the region in the "right" of the chosen edge is not
            % too small, add the edges of the current region to the edge_set field of
            % the second child of the current node
            if ~emptydx
                region_index1 = [region_index1 j];
                edge_set1 = [edge_set1 edges_to_add];
            end
        end
        
        % Add the two children to nodeset array
        nodeset = [nodeset child1 child2];
        nodeset = unique(nodeset);
        
        % Fill the children nodes fields
        nodes(child1).name = (nodes(ii).name-1)*2+2;
        nodes(child2).name = (nodes(ii).name-1)*2+3;
        nodes(child1).edge_index = edge_index;
        nodes(child2).edge_index = edge_index;
        nodes(child1).region_index = region_index1;
        nodes(child2).region_index = region_index2;
        nodes(child1).edge_set = unique(edge_set1);
        nodes(child2).edge_set = unique(edge_set2);
        nodes(child1).depth = nodes(ii).depth+1;
        nodes(child2).depth = nodes(ii).depth+1;
        % Remove chosen edge from edge_set
        nodes(child1).edge_set(nodes(child1).edge_set == edge_index) = [];
        nodes(child2).edge_set(nodes(child2).edge_set == edge_index) = [];
        
        if isempty(nodes(child1).edge_set) && numel(nodes(child1).region_index) > 1
            % If there are no more edges to choose for the first child and the number
            % of regions still available is greater than 1, the child is a leaf and
            % the first region is considered.
            nodes(child1).leaf = 1;
            nodes(child1).children = [];
            nodes(child1).region_index = nodes(child1).region_index(1);
            depths = [depths nodes(child1).depth];
        elseif numel(nodes(child1).region_index) == 1
            % If only one region is available, the child is a leaf and the region is
            % chosen
            nodes(child1).leaf = 1;
            nodes(child1).children = [];
            depths = [depths nodes(child1).depth];
        elseif size(unique([regions(nodes(child1).region_index).Ifunctions]','rows'),1) <= 1
            % If all the available regions have the same function, the child is a leaf
            % and the first region is chosen.
            nodes(child1).leaf = 1;
            nodes(child1).children = [];
            try
                nodes(child1).region_index = nodes(child1).region_index(1);
            catch
                nodes(child1).region_index = nodes(ii).region_index(1);
            end
            depths = [depths nodes(child1).depth];
        else
            % In the other cases, the child is not a leaf
            nodes(child1).leaf = 0;
            nodes(child1).children = [nodeset(end)+[1 2]];
            % ndeset is updated with the index of the child's children
            nodeset = [nodeset nodes(child1).children];
            nodeset = unique(nodeset);
        end
        
        
        % Same comments as for child 1
        if isempty(nodes(child2).edge_set) && numel(nodes(child2).region_index) > 1
            nodes(child2).leaf = 1;
            nodes(child2).children = [];
            nodes(child2).region_index = nodes(child2).region_index(2);
            depths = [depths nodes(child2).depth];
        elseif numel(nodes(child2).region_index) == 1
            nodes(child2).leaf = 1;
            nodes(child2).children = [];
            depths = [depths nodes(child2).depth];
        elseif size(unique([regions(nodes(child2).region_index).Ifunctions]','rows'),1) <= 1
            nodes(child2).leaf = 1;
            nodes(child2).children = [];
            try
                nodes(child2).region_index = nodes(child2).region_index(1);
            catch
                nodes(child2).region_index = nodes(ii).region_index(2);
            end
            depths = [depths nodes(child2).depth];
        else
            nodes(child2).leaf = 0;
            nodes(child2).children = [nodeset(end)+[1 2]];
            nodeset = [nodeset nodes(child2).children];
            nodeset = unique(nodeset);
        end
    end
    
    % If there are no more nodes to create, exit the loop
    if inodeset == numel(nodeset)
        proceed = 0;
    end
    inodeset = inodeset + 1;
    
end

close(h);

% Update tree fields
tree.numberOfNodes = numel(nodes);
tree.minDepth = min(depths);
tree.meanDepth = mean(depths);
tree.maxDepth = max(depths);
tree.nodes = nodes;
tree.balance = var(depths);

object(1).tree = tree;

if numel(object) > 1
    disp('Binary search tree computed for first element of pwag array:');
else
    disp('Binary search tree computed:');
end
disp(['   number of nodes: ',num2str(tree.numberOfNodes)]);
disp(['   minimum depth: ',num2str(tree.minDepth)]);
disp(['   mean depth: ',num2str(tree.meanDepth)]);
disp(['   maximum depth: ',num2str(tree.maxDepth)]);
disp(['   balance: ',num2str(tree.balance)]);
end