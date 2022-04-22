function regions = getRegions(object,varargin)
% getRegions  Gets information about the requested polytope(s).
%
% REGIONS = getNumberOfRegions(OBJ)
% Returns an array of structures (REGIONS) with as many elements as the
% number of polytopes composing the partition. Each elements of array 
% REGIONS is a structure with fields H, K, F and G. H and K define the 
% edges of the polytope, in the form:
%
%       H x <= K
%
% F and G define the affine function which lies over the polytope, in
% the form:
%
%       u = Fx+G
%
% OBJ is the pwagFunction object.
%
% REGIONS = getNumberOfRegions(OBJ,REG)
% Returns an array of structures (REGIONS) with as many elements as the
% number of elements in vector REG. REG is a vector containing the indices
% of the polytopes to get.

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



if nargin == 1
    ireg = 1:object.nr;
elseif nargin == 2
    ireg = varargin{1};
else
    error('To many input arguments');
end

nireg = numel(ireg);
regions = struct('H',cell(nireg,1),'K',cell(nireg,1),...
    'F',cell(nireg,1),'G',cell(nireg,1));

for i = 1:nireg
    
    % Indices of edges
    idx = object.regions(ireg(i)).Iedges(:,1);
    
    % Edges H
    Htmp = object.edges.H(idx,:);
    
    % Sign of the edges
    sign = object.regions(ireg(i)).Iedges(:,2);
    
    % Multiply by the sign
    Htmp = Htmp.*repmat(sign,1,object.nx);
    
    % Update field
    regions(i).H = Htmp;
    
    % Edges K
    Ktmp = object.edges.K(idx);
    
    % Multiply by the sign
    Ktmp = Ktmp.*sign;
    
    % Update field
    regions(i).K = Ktmp;
    
    % Indices of functions
    idx = object.regions(ireg(i)).Ifunctions;
    
    % Function F
    Ftmp = object.functions.F(idx,:);
    
    % Update field
    regions(i).F = Ftmp;
    
    % Function G
    Gtmp = object.functions.G(idx);
    
    % Update field
    regions(i).G = Gtmp;
    
end
end