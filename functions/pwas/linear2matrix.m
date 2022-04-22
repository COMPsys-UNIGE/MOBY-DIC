function matrixIdx = linear2matrix(linIdx,siz)
% linear2matrix   Compute the matrix indices starting from a linear index
%
% MIDX = linear2matrix(LIDX,SIZ)
% SIZ is an array with n elements specifying the size of the n-dimensional 
% matrix, i.e., for a 3 x 2 x 4 matrix, SIZ = [3 2 4]. LIDX is the linear
% index obtained by scanning the matrix by columns. MIDX is an array with n 
% elements with the indices for any matrix dimension.
%
% Example:
% Consider a 3 x 2 matrix.
% linear2matrix(1,[3 2]) = [1; 1]
% linear2matrix(2,[3 2]) = [2; 1]
% linear2matrix(3,[3 2]) = [3; 1]
% linear2matrix(4,[3 2]) = [1; 2]
% linear2matrix(5,[3 2]) = [2; 2]
% linear2matrix(6,[3 2]) = [3; 2]

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

% Number of dimension of the matrix
n = numel(siz);

if linIdx < 1 || linIdx > prod(siz)
    error('Wrong linear index.');
end

matrixIdx = zeros(n,1);
num = linIdx-1;
den = prod(siz(1:end-1));
for k = n:-1:2
    matrixIdx(k) = floor(num/den);
    num = num-matrixIdx(k)*den;
    den = den/siz(k-1);
end
matrixIdx(1) = floor(num/den);

matrixIdx = matrixIdx+1;
