function disp(object)
% disp   Displays some information about the embeddedSystem object
%
% disp(OBJ)
% OBJ is the MPCctrl object.

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
% Boston, MA  02111-1307  US

if object.hasObserver() && object.hasController()
    disp('Embedded system object containing controller');
    if isa(object.dynSys,'pwaSys')
        disp('and observer for a PWA system with: ');
    else
        disp('and observer for a LTI system with: ');
    end
    disp([' - ',num2str(object.nx),' state variables']);
    disp([' - ',num2str(object.nu),' input variables']);
    disp([' - ',num2str(object.np),' parameters']);
    disp([' - ',num2str(object.nd),' unmeasurable inputs']);
    if isa(object.dynSys,'pwaSys')
        disp([' - ',num2str(object.dynSys.getNumberOfDynamics),' dynamics']);
    end
    disp(' ')
    disp(['Controller sampling time: ',num2str(object.getController.getSamplingTime),' seconds']);
    disp(['Observer sampling time: ',num2str(object.getObserver.getSamplingTime),' seconds']);
    disp(' ')
    disp('Variation range of the variables are:')
    for i=1:object.nx
        disp([num2str(object.range.xmin(i)),' <= x',num2str(i),' <= ',num2str(object.range.xmax(i))]);
    end
    for i=1:object.np
        disp([num2str(object.range.pmin(i)),' <= p',num2str(i),' <= ',num2str(object.range.pmax(i))]);
    end
    for i=1:object.nd
        disp([num2str(object.range.dmin(i)),' <= d',num2str(i),' <= ',num2str(object.range.dmax(i))]);
    end
    for i=1:object.ny
        disp([num2str(object.range.ymin(i)),' <= y',num2str(i),' <= ',num2str(object.range.ymax(i))]);
    end
    for i=1:numel(object.range.xrefmin)
        disp([num2str(object.range.xrefmin(i)),' <= xref',num2str(i),' <= ',num2str(object.range.xrefmax(i))]);
    end
    for i=1:object.nu
        disp([num2str(object.range.umin(i)),' <= u',num2str(i),' <= ',num2str(object.range.umax(i))]);
    end
elseif object.hasObserver()
    disp('Embedded system object containing one observer.');
    object.getObserver().disp();
else
    disp('Embedded system object containing one controller.');
    object.getController().disp();
end
