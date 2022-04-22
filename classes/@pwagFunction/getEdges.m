function [H, K] = getEdges(object,varargin)
% getEdges  Gets the coefficients of the edges of a region
%
% An edge is defined as H x <= K.
%
% [H, K] = getEdges(OBJ)
% Returns all the edges of the polytopic partition. The edges are returned
% as a "big" matrix H of dimensions N_EDGES x NX and a "big"
% column vector K of N_EDGES elements. OBJ is the pwagFunction object.
% 
% [H, K] = getEdges(OBJ,REG)
% Returns the edges of the polytope whose index is REG.

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
    H = object.edges.H;
    K = object.edges.K;
elseif nargin == 2
    reg = varargin{1};
    
    if numel(reg) ~= 1
        error('REG must be a scalar!')
    end
    if reg < 0 || reg > object.nr
        error(['REG must be a scalar between 1 and ',num2str(object.nr)])
    end
    
    Iedges = object.regions(reg).Iedges;
    
    index = Iedges(:,1);
    sign = Iedges(:,2) == -1;
    
    H = object.edges.H(index,:);
    K = object.edges.K(index);
    H(sign,:) = -H(sign,:);
    K(sign) = -K(sign);
else
    error('Wrong input arguments');
end

end