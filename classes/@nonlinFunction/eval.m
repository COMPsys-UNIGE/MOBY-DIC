function u = eval(object,x,idx)
% eval   Evaluates the non linear function
%
% Y = eval(OBJ,X)
% Evaluates all NY components of the non linear function at the points
% specified by matrix X. X must be a NX x NPOINTS matrix. Y is a
% NY x NPOINTS matrix. NX and NY are the domain and codomain dimensions,
% respectively. OBJ is the nonlinFunction object.
%
% Y = eval(OBJ,X,IDX)
% Evaluates only the function components defined by IDX. IDX must be a vector
% of indices in the range 1, NY.

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
    error('IDX must be an integer number');
end
if any(idx < 1) || any(idx > object.ny)
    error(['IDX must be an integer number between 1 and ',num2str(object.ny)]);
end

% Check on x
if size(x,1) ~= object.nx
    error(['Matrix x must be a matrix with ',num2str(object.nx),' rows']);
end

% Get domain
Hd = object.domain.Hd;
Kd = object.domain.Kd;

% Number of points
npoints = size(x,2);

% Evaluate function
u = object.fhandle(x);
% Extract only requested components
u = u(idx,:);

% Tolerance
pars = getMOBYDICpars();
tol = pars.roundtol;

% Check if the point is outside the domain
ind = any(Hd*x > repmat(Kd,1,npoints)+tol);
u(:,ind) = NaN;

if any(isnan(u))
    warning('Function evaluation returned NaN')
end



end