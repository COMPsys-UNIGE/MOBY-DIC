function [r,xc] = chebyradius(H,K)
% chebyradius   Computes the Chebyshev radius and the Chebyshev center of a polytope
%
% [R,X] = chebyradius(H,K)
% Computes the Chebyshev radius R and the Chebyshev center X of a polytope
% defined by matrices H and K: H x <= K.
%
% MPT3 solver is used.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Tomaso Poggi (tpoggi@essbilbao.org)
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

P = Polyhedron(H,K);
s = P.chebyCenter;
r = s.r;
xc = s.x;
