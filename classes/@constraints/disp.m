function disp(object)
% disp   Displays some information about the constraints object
%
% disp(OBJ)
% OBJ is the constraints object.

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


if object.nx == 0
    disp('Empty constraints object');
else
    disp('Linear inequality constraints involving: ')
    if object.nx > 0
        disp(['- ',num2str(object.nx),' states'])
    end
    if object.nu > 0
        disp(['- ',num2str(object.nu),' inputs'])
    end
    if object.ny > 0
        disp(['- ',num2str(object.ny),' outputs'])
    end
    if object.np > 0
        disp(['- ',num2str(object.np),' parameters'])
    end
    if object.nd > 0
        disp(['- ',num2str(object.nd),' unmeasurable inputs'])
    end
    if object.nref > 0
        disp(['- ',num2str(object.nref),' references'])
    end
    disp(' ')
    disp(['Time horizon: ',num2str(object.N)]);
end
