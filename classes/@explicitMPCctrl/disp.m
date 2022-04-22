function disp(object)
% disp   Displays some information about the explicitMPCctrl object
%
% disp(OBJ)
% OBJ is the explicitMPCctrl object.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Mateo Lodi (matteo.lodi@edu.unige.it)
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
    disp('Empty explicitMPCctrl object');
else
    disp('Explicit MPC controller for a system with: ')
    disp([' - ',num2str(object.nx),' state variables']);
    disp([' - ',num2str(object.nu),' input variables']);
    disp([' - ',num2str(object.ny),' output variables']);
    disp([' - ',num2str(object.np),' parameters']);
    disp([' - ',num2str(object.nd),' unmeasurable inputs']);
    disp(' ')
    disp(['Sampling time: ',num2str(object.Ts),' seconds']);
    disp(' ')
    disp('The controller is a PWA function defined over a ')
    disp([num2str(object.getNumberOfDimensions()),'-dimensional domain partitioned into ',...
        num2str(object.getNumberOfRegions()),' polytopes']);
    disp(' ')
    
    if ~object.tracking
        if strcmpi(object.ctrlvar,'state')
            disp(['The controller is designed for the regulation of the system',...
                ' states to a fixed reference ',num2str(object.ref(:)')]);
        else
            disp(['The controller is designed for the regulation of the system',...
                ' outputs to a fixed reference ',num2str(object.ref(:)')]);
        end
    else
        if strcmpi(object.ctrlvar,'state')
            statenames = '';
            for i = 1:numel(object.trackvar)-1
                statenames = [statenames, '''',object.xnames{object.trackvar(i)},''', '];
            end
            statenames = [statenames, '''',object.xnames{object.trackvar(end)},''''];
            
            disp('The controller is designed for tracking')
            disp(['state variable(s) ',statenames]);
        else
            outputnames = '';
            for i = 1:numel(object.trackvar)-1
                outputnames = [outputnames, '''',object.ynames{object.trackvar(i)},''', '];
            end
            outputnames = [outputnames, '''',object.ynames{object.trackvar(end)},''''];
            
            disp('The controller is designed for tracking')
            disp(['output variable(s) ',outputnames]);
        end
    end
    disp(' ');
    
    tree = object.fun.getTree();
    if ~isempty(tree)
        disp('The tree is computed:');
        disp([' - number of nodes: ',num2str(tree.numberOfNodes)]);
        disp([' - minimum depth: ',num2str(tree.minDepth)]);
        disp([' - mean depth: ',num2str(tree.meanDepth)]);
        disp([' - maximum depth: ',num2str(tree.maxDepth)]);
        disp([' - balance: ',num2str(tree.balance)]);
    end
end