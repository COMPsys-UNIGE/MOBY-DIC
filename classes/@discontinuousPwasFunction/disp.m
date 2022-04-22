%% disp  displays some information about the discontinuous pwasFunction object

% disp(OBJ)
% OBJ is the discontinuous pwasFunction object.

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

function disp(object)

if object.nx == 0
    disp('Empty discontinuous pwasFunction object');
    disp(' ')
else
    disp('discontinuous pwasFunction object');
    disp(' ')
    disp([' - Domain dimensions: ',num2str(object.nx)]);

        
    xstr = [];
    for i = 1:object.nx-1
        xstr = [xstr object.xnames{i} ', '];
    end
    xstr = [xstr object.xnames{end}];
    ystr = [];
    for i = 1:object.ny-1
        ystr = [ystr object.ynames{i} ', '];
    end
    ystr = [ystr object.ynames{end}];
    disp(' ')
    disp([' - Input names: ',xstr])
    disp([' - Output names: ',ystr])
    disp(' ')
    
end
end