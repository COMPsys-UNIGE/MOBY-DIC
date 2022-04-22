function generateC(object,varargin)
% generateC
% Generates the C files describing the embedded system (observer + controller)
%
% generateC(OBJ)
% Generates the C files for the implementation on a microcontroller of the 
% embedded system with the default options.
% Several C files are generated, the main file is main.c.
% A report is also generated showing the main circuit features (latency,
% area occupation, etc.). A floating point (single) data representation 
% is used.
%
% generateC(OBJ,OPTS)
% Generates the C files for the implementation on a microcontroller of the 
% embedded system with custom options.
% OPTS is a structure with the following fields:
%  -  architecture: can be either 'fast' (default) or 'small'. The first
%                   choice produces a circuit with lower latency and bigger
%                   size, the second one viceversa.
%  -  inputResolution: number of bits used to represent the function
%                      inputs. In general, it corresponds to the resolution
%                      of the analog-to-digital converters (ADC), if any.
%                      Default value: 12.
%  -  inputRepresentation: can be either 'signed' or 'unsigned'. This also
%                          depends on the ADCs. Default value: 'signed'.
%  -  inputRange: structure with fields min and max, indicating the minimum
%                 and maximum value of the circuit inputs. min and max must
%                 be integer values representable with the number of bits
%                 chosen in inputResolution. For example if inputResolution
%                 = 8 and inputRepresentation = 'unsigned', min and max
%                 must be integer numbers comprised betwen 0 and 255. In
%                 this case, if 0 and 255 are chosen as inputRange.min and
%                 inputRange.max, the maximum representation range is
%                 exploited. Default: maximum representable range.
%  -  coeffResolution: number of bits used to represent every signal inside
%                      the algorithm used to solve the optimization
%                      problem. This value is in general bounded by the 
%                      resolution of the multipliers available on the FPGA.
%                      Default value: inputResolution.
%  -  coeffIntResolution: number of bits used to represent the integer part
%                      of the signal with coeffResolution bits. To be used
%                      only in case of implicit MPC controller.
%  -  outputResolution: number of bits used to represent the function
%                       outputs. In general, it corresponds to the 
%                       resolution of the digital-to-analog converters
%                       (DAC), if any. Default value: coeffResolution.
%  -  outputRepresentation: can be either 'signed' or 'unsigned'. This also
%                           depends on the DACs. Default value: 'signed'.
%  -  outputRange: structure with fields min and max, indicating the
%                  minimum and maximum value of the circuit outputs.
%                  Default: maximum representable range.
%  -  defaultOutput: control value that the controller must produce when
%                    the reset signal is low.
%  -  ADMMparameters: struct containing 'maxIter' and 'regPar', which are
%                     the number of iterations of the ADMM algorithm and a
%                     regularization parameter, respectively. Default
%                     values: 40, 2. To be used only in case of implicit
%                     MPC controller
%  -  folder: folder where C files are created. Default value: progressive
%             folder.
% A report is also generated showing the main circuit features (latency,
% area occupation, etc.).
%
% PERF = generateVHDL(OBJ,...)
% If an output argument is provided, the report is not automatically opened
% and a structure PERF containing some information about the circuit 
% implementation is returned. PERF is a structure with the following fields:
% TO DO
%
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
    circuit_parameters = [];
elseif nargin == 2
    circuit_parameters = varargin{1};
else
    error('Wrong input arguments');
end

if ~(object.hasObserver || object.hasController)
    error('Cannot create VHDL of an embedded system without a observer or a controller');
end

circuit_parameters = object.embeddedSystemCset(circuit_parameters);
   
if object.hasObserver && object.hasController
    % Generate embedded system with controller and observer     
    if isa(object.getController, 'explicitMPCctrl')
        object.explicitCtrlGenerateC(circuit_parameters);
    else
        object.implicitCtrlGenerateC(circuit_parameters);
    end    
    
elseif object.hasObserver
    % Generate only observer
    warning('Main interface generated through observer generateC method');
    obs = object.getObserver();
    circuit_parameters.generateInterface = true;
    obs.generateC(circuit_parameters);
    
elseif object.hasController
    % Generate only controller
    warning('Main interface generated through controller generateC method');
    ctrl = object.getController();
    circuit_parameters.generateInterface = true;
    ctrl.generateC(circuit_parameters);
end


