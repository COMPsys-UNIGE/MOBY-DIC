function disp(object)
% disp   Displays some information about the ApproxMPCctrl object
%
% disp(OBJ)
% OBJ is the ApproxMPCctrl object.

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
% Boston, MA  02111-1307  US

if object.nu == 0
    disp('Empty MPCctrl object');
else
    disp('Explicit approximate MPC controller for a system with: ')
    disp([' - ',num2str(object.nx),' state variables']);
    disp([' - ',num2str(object.nu),' input variables']);
    disp([' - ',num2str(object.np),' parameters']);
    disp([' - ',num2str(object.nd),' unmeasurable inputs']);
    disp(' ')
    disp(['Sampling time: ',num2str(object.Ts),' seconds']);
    disp(' ')
    disp('The controller is a PWAS function defined over a ')
    disp([num2str(object.getNumberOfDimensions()),'-dimensional domain partitioned into ',...
        num2str(object.getNumberOfRegions()),' simplices']);
    disp(' ')
    
    if ~object.tracking
        disp('The controller is designed for the regulation')
        disp('of the system state to a fixed reference')
    else
        statenames = '';
        for i = 1:numel(object.trackvar)-1
            statenames = [statenames, '''',object.xnames{object.trackvar(i)},''', '];
        end
        statenames = [statenames, '''',object.xnames{object.trackvar(end)},''''];
        
        disp('The controller is designed for tracking')
        disp(['state variable(s) ',statenames]);
    end
    disp(' ');
end