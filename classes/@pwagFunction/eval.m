function u = eval(object,x,idx)
% eval   Evaluates the PWAG function
%
% Y = eval(OBJ,X)
% Evaluates all NY components of the PWAG function at the points
% specified by matrix X. X must be a NX x NPOINTS matrix. Y is a 
% NY x NPOINTS matrix. NX and NY are the domain and codomain dimensions, 
% respectively. OBJ is the pwagFunction object.
%
% Y = eval(OBJ,X,IDX)
% Evaluates only the function components defined by IDX. IDX must be a vector
% of indices in the range 1, NY.
%
% Mex files are used to speedup the computation.

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

if ~exist('idx','var')
    idx = 1:object.ny;
end

if floor(idx) ~= idx
    error('idx must be an integer number');
end
if any(idx < 1) || any(idx > object.ny)
    error(['idx must be an integer number between 1 and ',num2str(object.ny)]);
end

ny = numel(idx);

% Check on x
if size(x,1) ~= object.nx
    error(['Matrix x must be of a matrix with ',num2str(object.nx),' rows']);
end

% Number of points x
npoints = size(x,2);

usemex = 1;

if ~usemex
    
    % Find region
    reg = object.findRegion(x);
    
    u = zeros(ny,npoints);
    for i = 1:npoints
        xcur = x(:,i);
        if reg(i) == 0
            u(:,i) = NaN;
        else
            [F, G] = object.getFunctions(reg(i));
            u(:,i) = F(idx,:)*xcur+G(idx);
        end
    end
    
else
    
    % Get matrices defining the edges
    [H, K] = object.getEdges();
    
    % Number of regions
    nr = object.nr;
    
    % Number of edges
    nH = size(H,1);
    
    % Initialize arrays
    nindexes = zeros(nr,1);
    indexes = zeros(nr*nH,1);
    ifunctions = zeros(nr*ny,1);
    signs = zeros(nr*nH,1);
    offset = 0;
    
    % Fill arrays
    % Loop on regions
    for i = 1:nr
        tmp = object.regions(i).Iedges;
        % Number of edges of i-th region
        ntmp = size(tmp,1);
        % Number of edges of each region
        nindexes(i) = ntmp;
        % Indexes of the edges of each region
        indexes(offset+[1:ntmp]) = tmp(:,1);
        % Signs of the edges of each region
        signs(offset+[1:ntmp]) = tmp(:,2);
        offset = offset+ntmp;
    end
    
    % Delete exceeding entries
    indexes(offset+1:end) = [];
    signs(offset+1:end) = [];
    
    % Take the matrices F and G rearranged as  column vectors
    F = object.functions.F;
    G = object.functions.G;
    
    % Fill an array ifunctions with the indexes of the functions
    % to compute for each region
    for i = 1:nr
        tmp = object.regions(i).Ifunctions(idx);
        ifunctions(ny*(i-1)+1:ny*(i-1)+ny) = tmp;
    end
    
    [Hd, Kd] = object.getDomain();
    
    % Call mex function
    [u, regmex] = pwagFunction_evalmex(ny,nr,H,K,F,G,nindexes,indexes,ifunctions,signs,Hd,Kd,x');    
    % The mex function returns 1.0e99 if the point is out of the
    % domain. This value is transformed in NaN
    u(u > 1.0e90) = NaN;
       
end