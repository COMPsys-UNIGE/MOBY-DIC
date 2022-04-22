function [Hm,Km,empty] = minimalRepresentation(H,K)
% minimalRepresentation   Finds the minimal representation of a polytope
%
% [Hm,Km,EMPTY] = minimalRepresentation(H,K)
% Finds the minimal representation Hm x <= Km of a polytope defined by
% H x <= K. EMPTY is set to true if the polytope is empty.
%
% MPT3 solver is used.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
%
% Copyright (C) 2015 University of Genoa, Italy.

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

% Create polyhedron
P = Polyhedron(H,K);

% Minimal representation
P.minHRep();

if P.isFullDim
    
    % Go back to H and K
    HH = P.H;
    
    Hm = HH(:,1:end-1);
    Km = HH(:,end);
    
    empty = 0;
    
else
    
    Hm = [];
    Km = [];
    empty = 1;
end

