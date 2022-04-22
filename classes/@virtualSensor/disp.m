function disp(object)
% disp   Displays some information about the virtualSensor object
%
% disp(OBJ)
% OBJ is the virtualSensor object.

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

if object.nz == 0
    disp('Empty virtualSensor object');
    disp(' ')
else
    disp('virtualSensor object');
    disp(' ')
    disp([' - Number of inputs: ',num2str(object.nu)]);
    disp([' - Number of measurable outputs: ',num2str(object.ny)]);
    disp([' - Number of unmeasurable outputs: ',num2str(object.nz)]);
    disp([' - Input time window: ',num2str(object.mu)]);
    disp([' - Output time window: ',num2str(object.my(:)')]);
    disp([' - Autoregressive time window: ',num2str(object.mz)]);
    
    ustr = [];
    if object.nu > 0
        for i = 1:object.nu-1
            ustr = [ustr object.unames{i} ', '];
        end
        ustr = [ustr object.unames{end}];
    end
    ystr = [];
    if object.ny > 0
        for i = 1:object.ny-1
            ystr = [ystr object.ynames{i} ', '];
        end
        ystr = [ystr object.ynames{end}];
    end
    zstr = [];
    for i = 1:object.nz-1
        zstr = [zstr object.znames{i} ', '];
    end
    zstr = [zstr object.znames{end}];
    disp(' ')
    disp([' - Input names: ',ustr])
    disp([' - Measurable output names: ',ystr])
    disp([' - Unmeasurable output names: ',zstr])
    disp(' ')
    if object.current
        disp('u and y at current time instant are used to estimate z')
    else
        disp('u and y at current time instant are not used to estimate z')
    end
    disp(' ')
    if object.identified
        disp('The virtual sensor has been identified.');
        disp(['Number of subdivisions per dimension: ',num2str(object.np)]);
    else
        disp('The virtual sensor is not identified.');
    end
    disp(' ')
    
end
