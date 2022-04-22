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
% generateSimulinkModel(OBJ,OPTS,simulateVHDL)
% generateSimulinkModel(OBJ,OPTS,simulateVHDL,circuit_parameters)
% OPTS is a structure with the following fields:
%
%   - folder: the folder in which the model is created
%
% simulateVHDL : if true the VHDL files describing the digital circuit
%                 implementing the  controller and/or observer are simulated. 
%                 This requires Xilinx System Generator. 
%                 If false, the simulation employs simply the
%                 controller and/or observer object. Default value false.
%
% circuit_parameters: is a structure describe the circuit parameters 
%                       with the following fields:
%
%       * inputResolution: number of bits used to code inputs (u,p,y) (default
%         inputResolution = 12)
%       * inputRange: variation range of the input. It is a struct with field
%               min and max (default full range)
%       * inputRepresentation: representation of the input (signed,unsigned)
%                       (default signed)
%
%       * coeffResolution: number of bits used to code the coefficients in memory
%                   (default coeffResolution = 12)
%
%       * outputResolution: number of bits used to code outputs (x,d) (default
%         outputResolution = 12)
%       * outputRange: variation range of the output. It is a struct with field
%               min and max (default full range)
%       * outputRepresentation: representation of the output (signed,unsigned)
%                       (default signed)
%
%       * frequency: indicates the frequency (in Hz) at which the circuit will
%              work. Default frequency = 50000000 (50 MHz).
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


if nargin == 1
    folder = '';
    simulateVHDL = false;
    circuit_parameters = [];
elseif nargin == 2
    folder = varargin{1};
    simulateVHDL = false;
    circuit_parameters = [];
elseif nargin == 3
    folder = varargin{1};
    simulateVHDL = varargin{2};
    circuit_parameters = [];
elseif nargin == 4
    folder = varargin{1};
    simulateVHDL = varargin{2};
    circuit_parameters = varargin{3};
else
    error('Wrong number of input arguments');
end


if ~ischar(folder)
    error('OPTS.folder must be a string');
end

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

circuit_parameters.folder = [folder,'/VHDLfile'];

options.simulinkFolder = getsimulinkpath();
options.folder = folder;
options.simulateVHDL = simulateVHDL;
options.circuit_parameters = circuit_parameters;

% Simulation of continuous-time system
if object.isContinuousTime
    
    % Open-loop system
    if ~object.hasController()
        
        if ~object.hasObserver()
            
            % Without observer
            ct_simulinkOpenLoopSim(object,options);
            
        else
            
            % With observer
            ct_simulinkOpenLoopObsSim(object,options);
        end
        
        % Closed-loop system
    else
        
        % Without observer
        if ~object.hasObserver()
            
            ct_simulinkClosedLoopSim(object,options);
            
            % With observer
        else
            
            ct_simulinkClosedLoopObsSim(object,options);
            
        end
        
    end
    
    % Simulation of discrete-time system
else
    
    % Open-loop system
    if ~object.hasController()
        
        if ~object.hasObserver()
            
            % Without observer
            dt_simulinkOpenLoopSim(object,options);
            
        else
            
            % With observer
            dt_simulinkOpenLoopObsSim(object,options);
        end
        
        % Closed-loop system
    else
        
        % Without observer
        if ~object.hasObserver()
            
            dt_simulinkClosedLoopSim(object,options);
            
            % With observer
        else
            
            dt_simulinkClosedLoopObsSim(object,options);
            
        end
        
    end
    
end

end