function reg = findRegion(object,x)
% findRegion   Finds the polytope containing a given point
%
% REG = findRegion(OBJ,X)
% Finds the indices REG of the polytopes containing points X. X must be a
% matrix of dimension NX x N_POINTX, being NX the number of domain
% dimension. REG is a vector with N_POINTS elements. OBJ is the
% pwagFunction object.
%
% Mex files are used to speedup the computation

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

if size(x,1) ~= object.nx
    error(['Matrix x must be of a matrix with ',num2str(object.nx),' rows']);
end

npoints = size(x,2);    % Number of points
nr = object.nr;         % Number of regions

% Matrices defining the edges of the regions
H = object.edges.H;
K = object.edges.K;

% Domain
Hd = object.domain.Hd;
Kd = object.domain.Kd;

usemex = 1;

if ~usemex
    
    % Initialize output
    reg = zeros(npoints,1);
    
    % Tolerance used to prevent that points in the edges of a region are
    % considered outside it
    % Retrieve rounding tolerance
    pars = getMOBYDICpars();
    tol = pars.roundtol;
    
    for ii = 1:npoints
        xx = x(:,ii);
        
        inDomain = Hd*xx <= Kd+tol;
        
        % Search the region
        reg(ii) = 0;
        
        if inDomain
            for i = 1:nr
                Iedges = object.regions(i).Iedges;
                
                index = Iedges(:,1);
                sign = Iedges(:,2) == -1;
                
                Hloc = H(index,:);
                Kloc = K(index);
                Hloc(sign,:)=-Hloc(sign,:);
                Kloc(sign)=-Kloc(sign);
                
                found = Hloc*xx-Kloc <= tol;
                if found
                    reg(ii) = i;
                    break
                end
            end
        end
    end
    
else
    
    % Number of edges
    nH = size(H,1);
    
    % Initialize arrays
    nindexes = zeros(nr,1);
    indexes = zeros(nr*nH,1);
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
    
    reg = pwagFunction_findRegionmex(nr,H,K,nindexes,indexes,signs,Hd,Kd,x);
    
end