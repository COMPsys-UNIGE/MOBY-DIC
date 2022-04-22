function object = computeTree(object,varargin)
% computeTree   Computes the binary search tree associated to the polytopic partition
%
% OBJ = computeTree(OBJ)
% Computes the bynary search tree with the default settings (see below).
% OBJ is the MPCctrl object.
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
% See also: MPCctrl/plotTree.

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

object.fun = object.fun.computeTree(mode,balance);