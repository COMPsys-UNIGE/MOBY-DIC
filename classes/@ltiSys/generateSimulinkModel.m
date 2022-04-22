function generateSimulinkModel(object,varargin)
% generateSimulinkModel    Generates a Simulink model for the system
%                          simulation
%
% generateSimulinkModel(OBJ)
% Generates a Simulink model for the simulation of the system.
% If an observer and/or a controller are associated to the system, they
% will be included in the model. OBJ is the ltiSys object.
%
% generateSimulinkModel(OBJ,OPTS)
% OPTS is a structure with the following fields:
%  - folder: the folder in which the model is created
%  - simVHDL: if true, the VHDL files implementing the controller and/or
%             observer associated to the system are generated and a
%             Simulink model allowing the hardware-in-the-loop simulation
%             (Simulink+FPGA) is created. This functionality requires 
%             Xilinx System Generator. If false, the simulation simply 
%             employs the controller and/or observer MATLAB objects. 
%             Default value false. 
%
% generateSimulinkModel(OBJ,OPTS,CIR_OPTS)
% If OPTS.simVHDL = 1, the options for the generation of the VHDL files can
% also be provided in structure CIR_OPTS. Refer to the documentation of 
% embeddedSystem.generateVHDL for a list of the structure fields.

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

% Check input arguments
if nargin == 1
    opts = '';
    cir_opts = '';
elseif nargin == 2
    opts = varargin{1};
    cir_opts = '';
elseif nargin == 3
    opts = varargin{1};
    cir_opts = varargin{2};
else
    error('Wrong number of input arguments');
end

if isempty(opts)
    folder = '';
elseif ~isfield(opts, 'folder')
    folder = '';
else
    folder = opts.folder;
end

if isempty(opts)
    simVHDL = 0;
elseif ~isfield(opts, 'simVHDL')
    simVHDL = 0;
else
    simVHDL = opts.simVHDL;
end

if ~isempty(cir_opts) && simVHDL == 0
    warning('Argument CIR_OPTS is ignored since OPTS.simVHDL = 0');
end

if ~ischar(folder)
    error('OPTS.folder must be a string');
end

% Create a folder
if ~isempty(folder)
    % Add \ at the end of the folder name
    if strcmp(folder(end),'\')
        folder(end) = '/';
    elseif ~strcmp(folder(end),'/')
        folder = [folder,'/'];
    end
else
    % Choose automatically a folder name
    created = 0;
    i = 1;
    while ~created
        testfolder = [pwd,'/model_',num2str(i),'/'];
        if ~exist(testfolder,'dir')
            folder = testfolder;
            created = 1;
        else
            i = i+1;
        end
    end
end

cir_opts.folder = [folder,'VHDLfile'];

options.simulinkFolder = getsimulinkpath();
options.folder = folder;
options.simulateVHDL = simVHDL;
options.circuit_parameters = cir_opts;

% Open-loop system
if ~object.hasController()

    if ~object.hasObserver()

        % Without observer
        simulinkOpenLoopSim(object,options);

    else

        % With observer
        simulinkOpenLoopObsSim(object,options);
    end

% Closed-loop system
else

    % Without observer
    if ~object.hasObserver()

        simulinkClosedLoopSim(object,options);

    % With observer
    else
    
        simulinkClosedLoopObsSim(object,options);

    end

end

end