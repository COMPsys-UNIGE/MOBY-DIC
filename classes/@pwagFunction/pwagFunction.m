classdef pwagFunction < MOBYFunction
    % pwagFunction   Piece-Wise Affine Generic function
    %
    % This object represents a piecewise affine function defined over a
    % generic polytopic domain partition. A PWAG function is defined as
    % follows:
    %
    %       u = F_i x + G_i,   if x belongs to P_i
    %
    % being P_i convex polytopes and F_i and G_i matrices of appropriate 
    % dimension.
    %
    % OBJ = pwagFunction()
    % Builds an empty pwagFunction object OBJ.
    %
    % OBJ = pwagFunction(REGIONS)
    % Builds a pwagFunction object OBJ by specifying all the regions
    % composing the polytopic partition and the affine function defined
    % over them. REGIONS is an array of structures (one for each polytope)
    % with the following fields: H, K, F, G. H and K define the edges of the
    % polytope, in the form:
    %
    %       H x <= K
    %
    % F and G define the affine function which lies over the polytope, in
    % the form:
    %
    %       u = Fx+G
    %
    % The function domain is automatically computed by joining all regions
    % through MPT3.
    %
    % OBJ = pwagFunction(REGIONS,DOMAIN)
    % Builds a pwagFunction object OBJ by specifying all the regions and
    % the domain the function is defined over.
    %
    % If the domain is hyper rectangular, DOMAIN can be a
    % structure with fields xmin and xmax, providing the domain boundaries:
    %
    %    xmin <= x <= xmax
    %
    % Otherwise DOMAIN must be a struct with fields Hd and Kd, which define
    % the polytopic domain as follows:
    %
    %    Hd x <= Kd
    %
    % pwagFunction methods:
    %   computeTree - computes the binary search tree associated to the polytopic partition
    %   disp - displays some information about the pwagFunction object.
    %   eval - evaluates the PWAG function.
    %   findRegion - finds the polytope containing a given point.
    %   generateC - Generates C files for the circuit implementation of the PWAG function on microcontroller
    %   generateVHDL - Generates VHDL files for the circuit implementation of the PWAG function on FPGA
    %   getEdges - gets the coefficients of the edges of a region.
    %   getFunctions - gets the coefficients of the affine function defined over a region.
    %   getNumberOfEdges - gets the number of different edges of the domain partition.
    %   getNumberOfRegions - gets the number of polytopes of the domain partition.
    %   getRegions - gets information about the requested polytope(s).
    %   getTree - gets the binary search tree associated to the domain partition
    %   getVertices - gets the vertices of the polytopes of the domain partition.
    %   hasTree - Returns true if the binary search tree for the pwagFunction has been computed
    %   plot - plots the PWAG function.
    %   plotPartition - plots the polytopic domain partition.
    %   plotTree - plots the binary search tree associated to the domain partition.
    %
    % The pwagFunction object is derived from MOBYFunction and inherits
    % all its methods.
    %
    % See also MOBYFunction, pwasFunction.
    
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
    
    
    % Properties
    
    properties (Access = private)
        nr = [];
        edges = struct('H',[],'K',[]);
        functions = struct('F',[],'G',[]);
        regions = struct('Iedges',[],'Ifunctions',[]);
        tree = [];
    end
    
    
    % Methods
    
    methods
        
        % Constructor
        function object = pwagFunction(varargin)
            
            % Empty object
            if nargin == 0
                object.nx = [];
                object.ny = [];
                object.domain.Hd = [];
                object.domain.Kd = [];
                object.xnames = [];
                object.ynames = [];
                object.nr = [];
                object.edges = struct('H',[],'K',[]);
                object.functions = struct('F',[],'G',[]);
                object.regions = struct('Iedges',[],'Ifunctions',[]);
                object.tree = [];
                
            elseif nargin == 1 || nargin == 2
                
                regs = varargin{1};
                
                if ~isfield(regs,'H') || ~isfield(regs,'K')...
                        || ~isfield(regs,'F') || ~isfield(regs,'G')
                    error('Structure "regions" must have fields ''H'', ''K'', ''F'' and ''G''');
                end
                
                if nargin == 2
                    
                    domain = varargin{2};
                    
                    if isfield(domain,'xmin')
                        if isfield(domain,'Hd') || isfield(domain,'Kd')
                            error('domain must be a struct with fields "xmin, xmax" or "Hd, Kd"');
                        end
                        if ~isfield(domain,'xmax')
                            error('domain must be a struct with fields xmin, xmax or Hd, Kd');
                        end
                    elseif isfield(domain,'Hd')
                        if isfield(domain,'xmin') || isfield(domain,'xmax')
                            error('domain must be a struct with fields "xmin, xmax" or "Hd, Kd"');
                        end
                        if ~isfield(domain,'Kd')
                            error('domain must be a struct with fields xmin, xmax or Hd, Kd');
                        end
                    end
                    

                    
                    % Number of domain dimensions
                    nx = size(regs(1).H,2);
                    
                    if isfield(domain,'xmin')
                        
                        % Transform into column vectors
                        domain.xmin = domain.xmin(:);
                        domain.xmax = domain.xmax(:);
                        
                        if numel(domain.xmax) ~= nx || numel(domain.xmin) ~= nx
                            error(['domain.xmin and domain.xmax must have ',num2str(nx),' elements']);
                        end
                        
                        if any(domain.xmin >= domain.xmax)
                            error('domain.xmin must be lower than domain.xmax');
                        end
                        
                        % Create matrices to express the domain
                        domain.Hd = [eye(nx);-eye(nx)];
                        domain.Kd = [domain.xmax;-domain.xmin];
                    end
                    
                    Hd = domain.Hd;
                    Kd = domain.Kd;
                    
                    % check if domain is bounded
                    pol = Polyhedron(Hd,Kd);
                    if ~pol.isBounded
                        error('Function domain must be bounded');
                    end
                    
                    if pol.isEmptySet
                        error('Domain is empty');
                    end
                    
                else
                    domain = [];
                    Hd = [];
                    Kd = [];
                end
                
                if size(Hd,1) ~= size(Kd,1)
                    error('domain.Hd and domain.Ks must have the same n umber of rows')
                end
                
                %number of domain dimensions
                nx = size(regs(1).F,2);              
                % Number of codomain dimensions
                ny = size(regs(1).F,1);
                % Number of regions
                nr = numel(regs);
                
                % Check domain edges (if provided)
                if ~isempty(domain)
                    if size(Hd,2) ~= nx
                        error(['domain.Hd must have ',num2str(nx),' columns']);
                    end
                    if size(Kd,2) ~= 1
                        error('domain.Kd must have 1 column');
                    end
                end
                
                
                xnames = cell(nx,1);
                ynames = cell(ny,1);
                for i = 1:nx
                    xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:ny
                    ynames{i} = ['y_',num2str(i)];
                end
                
                % Initialize matrices containing H, K, F and G
                HKregion = cell(nr,1);
                FGregion = cell(nr,1);
                
                % Initialize indices pointing to correct values of H
                % and K
                i1 = zeros(nr,1);
                i2 = zeros(nr,1);
                i2prev = 0;
                
                
                k = 1;
                
                for i = 1:nr
                    
                    % Extract structure fields
                    H = regs(i).H;
                    K = regs(i).K;
                    F = regs(i).F;
                    G = regs(i).G;
                    
                    % Check fields
                    if size(H,2) ~= nx
                        error(['Region ',num2str(i),': matrix H must have ',num2str(nx),' columns']);
                    end
                    if size(K,2) ~= 1
                        error(['Region ',num2str(i),': matrix K must have one column']);
                    end
                    if size(F,2) ~= nx
                        error(['Region ',num2str(i),': matrix F must have ',num2str(nx),' columns']);
                    end
                    if size(G,2) ~= 1
                        error(['Region ',num2str(i),': matrix G must have one column']);
                    end
                    if size(F,1) ~= ny
                        error(['Region ',num2str(i),': matrix F must have ',num2str(ny),' columns']);
                    end
                    if size(G,1) ~= ny
                        error(['Region ',num2str(i),': matrix G must have ',num2str(ny),' columns']);
                    end
                    if size(H,1) ~= size(K,1)
                        error(['Region ',num2str(i),': matrices H and K must have the same number of rows']);
                    end
                    
                    % Find minimal representation
                    [Hmin,Kmin,empty] = minimalRepresentation([H;Hd],[K;Kd]);
                    
                    % Update edges and indexes
                    % Regions outside the domain are removed
%                     if ~empty
                        % Number of edges of current polytope
                        nedges = size(Hmin,1);
                        
                        % Update edges, functions and indexes
                        HKregion{k} = [Hmin Kmin];
                        FGregion{k} = [F G];
                        i1(k) = i2prev+1;
                        i2(k) = i1(k)+nedges-1;
                        i2prev = i2(k);
                        k = k+1;
                        
%                     else
%                         disp(['WARNING: region ',num2str(i),' is ignored because it is empty or outside the domain!']);
%                     end
                end
                
                % If not provided, compute domain (through MPT3)
                if isempty(domain)
                    % Create a polyhedron for each region
                    for i = 1:nr
                        P(i) = Polyhedron(regs(i).H,regs(i).K); %#ok<*AGROW>
                    end
                    % Compute PolyUnion
                    U = PolyUnion(P);
                    % Compute convex hull
                    HU = U.convexHull;
                    
                    % Extract inequality matrices
                    HKd = HU.H;
                    Hd = HKd(:,1:end-1);
                    Kd = HKd(:,end);
                    
                    % Create domain structure
                    domain.Hd = Hd;
                    domain.Kd = Kd;
                end
                
                % New number of regions
                nr = k-1;
                HK = cell2mat(HKregion(1:nr));
                FG = cell2mat(FGregion(1:nr));
                
                % Remove unused entries
                i1(k:end) = [];
                i2(k:end) = [];
            else
                error('Wrong input arguments.');
            end
            
            MOBYDICpars = getMOBYDICpars;
            tol = MOBYDICpars.roundtol;
            
            % Round matrices to prevent numerical errors
            HKround = round(HK/tol)*tol;
            
            % Find equal or opposite rows in matrices H and K and
            % their correspondances
            
            % correspHK is a 3 columns array. First column is the
            % index of a row of matrix HK, second element is
            % another row with the same (or opposite) values. Third
            % row is 1 if the elements are the same or -1 if they
            % are opposite
            correspHK = zeros(size(HK,1),3);
            
            k = 1;  % Index for correspHK
            for i = 1:size(HK,1)
                if ~any(correspHK(:,1) == i)
                    
                    % Fill first row
                    correspHK(k,:) = [i i 1];
                    k = k+1;
                    
                    % Find equal rows
                    iplus = find(ismember(HKround,HKround(i,:),'rows'));
                    
                    % Find opposite rows
                    iminus = find(ismember(HKround,-HKround(i,:),'rows'));
                    iminus = [i; iminus];
                    
                    if numel(iplus)~=1
                        correspHK(k:k+numel(iplus)-2,1)=iplus(2:end);
                        correspHK(k:k+numel(iplus)-2,2)=iplus(1);
                        correspHK(k:k+numel(iplus)-2,3)=1;
                        k = k+numel(iplus)-1;
                    end
                    
                    if numel(iminus)~=1
                        correspHK(k:k+numel(iminus)-2,1)=iminus(2:end);
                        correspHK(k:k+numel(iminus)-2,2)=iminus(1);
                        correspHK(k:k+numel(iminus)-2,3)=-1;
                        k = k+numel(iminus)-1;
                    end
                end
            end
            
            % Create new reduced H and K matrices without redundant
            % rows
            HKred = zeros(size(HK));
            k = 1;  % Index for HKred
            
            % The elements of array indexRelation are such as:
            % HK(i,:) = HKred(indexRelation(i),:)
            indexRelation = zeros(size(HK,1),1);
            for i = 1:size(HK,1)
                j = correspHK(i,1);
                if j==correspHK(i,2)
                    HKred(k,:) = HK(j,:);
                    indexRelation(j) = k;
                    k = k+1;
                end
            end
            
            % Delete empty rows, if any
            HKred(k:end,:) = [];
            
            % Update correspondances based on the new H and K
            % matrices
            correspHK(:,2) = indexRelation(correspHK(:,2));
            
            % Find equal rows in matrices F and G and
            % their correspondances.
            FGround = round(FG/tol)*tol;
            
            % correspFG is a 2 columns array. First column is the
            % index of a row of matrix FG, second element is
            % another row with the same values.
            correspFG = zeros(size(FG,1),2);
            
            % Index for correspFG
            k = 1;
            for i = 1:size(FG,1)
                if ~any(correspFG(:,1) == i)
                    
                    % Find equal rows
                    iplus = find(ismember(FGround,FGround(i,:),'rows'));
                    
                    correspFG(k,:) = [i i];
                    k = k+1;
                    
                    if numel(iplus)~=1
                        correspFG(k:k+numel(iplus)-2,1) = iplus(2:end);
                        correspFG(k:k+numel(iplus)-2,2) = iplus(1);
                        k = k+numel(iplus)-1;
                    end
                end
            end
            
            % Create new reduced F and G matrices without
            % redundant rows
            FGred = zeros(size(FG));
            
            % Index for FGred
            k = 1;
            
            % The elements of array indexRelation are such as:
            % FG(i,:) = FGred(indexRelation(i),:)
            indexRelation = zeros(size(FG,1),1);
            for i = 1:size(FG,1)
                j = correspFG(i,1);
                if j == correspFG(i,2)
                    FGred(k,:) = FG(j,:);
                    indexRelation(j,:) = k;
                    k = k+1;
                end
            end
            
            % Delete empty rows, if any
            FGred(k:end,:) = [];
            
            % Update correspondances based on the new F and G
            % matrices
            correspFG(:,2) = indexRelation(correspFG(:,2));
            
            % Create regions structure
            regions = struct('Iedges',cell(nr,1),'Ifunctions',cell(nr,1));
            
            for i = 1:nr
                Iedges = i1(i):i2(i);
                Ifunctions = ny*(i-1)+1:ny*(i-1)+ny;
                
                Iedges = correspHK(ismember(correspHK(:,1),Iedges),2:3);
                
                [ia, ib] = intersect(correspFG(:,1),Ifunctions); %#ok<ASGLU>
                Ifunctions = correspFG(ib,2);
                
                regions(i).Iedges = Iedges;
                regions(i).Ifunctions = Ifunctions;
            end
            
            % Fill object fields
            
            object.nx = nx;
            
            object.ny = ny;
            
            object.nr = nr;
            
            object.edges.H = HKred(:,1:end-1);
            object.edges.K = HKred(:,end);
            
            object.functions.F = FGred(:,1:end-1);
            object.functions.G = FGred(:,end);
            
            object.regions = regions;
            
            object.xnames = xnames;
            object.ynames = ynames;
            
            object.domain = domain;
            
        end
        
        % Build a PWA function given the regions and the function coefficients
        
        
        % Get methods
        nreg = getNumberOfRegions(object);
        
        nedges = getNumberOfEdges(object);
        
        regions = getRegions(object,varargin)

        tree = getTree(object)
        
        [H, K] = getEdges(object,varargin);
        
        [F, G] = getFunctions(object,varargin);
        
        V = getVertices(object,varargin);
        
        % Other methods
        
        object = computeTree(object,varargin)

        reg = findRegion(object,x,modality);
        
        y = eval(object,x,varargin);
        
        plot(object,varargin);
        
        plotPartition(object);

        plotTree(object);  

        answ = hasTree(object);  
        
        disp(object);

        varargout = generateVHDL(object,varargin);

        varargin = generateC(object,varargin);
        
    end
end


