function varargout = generateVHDL(object,varargin)
% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the embedded system (observer + controller)
%
% generateVHDL(OBJ)
% Generates the VHDL files for the circuit implementation on FPGA of the 
% embedded system with the default options.
% Several VHDL files are generated, the top-level block is
% embeddedSystemInterface.vhd. A report is also generated showing the main
% circuit features (latency, area occupation, etc.). A fixed point data
% representation is used.
%
% generateVHDL(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation on FPGA of the 
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
%                      of the signal with coeffResolution bits.
%  -  outputResolution: number of bits used to represent the function
%                       outputs. In general, it corresponds to the 
%                       resolution of the digital-to-analog converters
%                       (DAC), if any. Default value: coeffResolution.
%  -  outputRepresentation: can be either 'signed' or 'unsigned'. This also
%                           depends on the DACs. Default value: 'signed'.
%  -  outputRange: structure with fields min and max, indicating the
%                  minimum and maximum value of the circuit outputs.
%                  Default: maximum representable range.
%  -  useADC: flag indicating if input values come from an ADC. If it is 1,
%             input values are inside the range indicated by InputRange, 
%             otherwise input values are inside the range indicated by 
%             inputRange.
%  -  useDAC: flag indicating if output values go to a DAC. If it is 1, 
%             output values are inside the range indicated by outputRange,
%             otherwise output values are inside the range indicated by
%             outputRange.
%  -  defaultOutput: control value that the controller must produce when
%                    the reset signal is low.
%  -  frequency: clock frequency (in Hz) of the FPGA. This information is 
%                used to calculate the circuit latency.  Default value:
%                50000000 (50 MHz).
%  -  fpgaBoard: FPGA board target for the synthesis of C++ code with HLS.
%  -  ADMMparameters: struct containing 'maxIter' and 'regPar', which are
%                     the number of iterations of the ADMM algorithm and a
%                     regularization parameter, respectively. Default
%                     values: 40, 2.
%  -  generateInterface: flag indicating if a circuit interface must be
%                        generated. If the controller is part of an 
%                        embedded system, the interface must not be
%                        generated. Default value: true.
%  -  folder: folder where VHDL files are created. Default value:
%             progressive folder.
% A report is also generated showing the main circuit features (latency,
% area occupation, etc.).
%
% PERF = generateVHDL(OBJ,...)
% If an output argument is provided, the report is not automatically opened 
% and a structure PERF containing some information about the circuit 
% implementation is returned. PERF is a structure with the following fields:
%  - range: structure with the following fields:
%            - xmin, xmax: actual ranges of the function inputs
%            - umin, umax: actual ranges of the function outputs
%            - xrefmin, xrefmax: actual ranges of the reference signals
%            - ymin, ymax: actual ranges of the system outputs
%            - pmin, pmax: actual ranges of the system parameters
%            - dmin, dmax: actual ranges of the system unmeasurable inputs
%  - circuit_paremeters: parameters of the generated circuit
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

%% Generate embedded system with controller and observer
if object.hasObserver && object.hasController
    
    % Generate the VHDL code directly, if the MPC controller is explicit
    if ~isa(object.getController, 'implicitMPCctrl')

        optO = explicitCtrlGenerateVHDL(object, circuit_parameters);
        
        if nargout > 1
        error('Wrong number of outputs')
        elseif nargout == 1
        varargout{1} = optO;
        end
     
    % Generate the VHDL code through Vitis HLS, if the MPC controller is implicit
    else
        
        optO = implicitCtrlGenerateVHDL(object, circuit_parameters);
        
        if nargout > 1
        error('Wrong number of outputs')
        elseif nargout == 1
        varargout{1} = optO;
        end
        
    end
    
elseif object.hasObserver
    % generate only observer
    warning('Main interface generated through observer generateVHDL method');
    obs = object.getObserver();
    circuit_parameters.generateInterface = true;
    obs.generateVHDL(circuit_parameters,es.range);
elseif object.hasController
    % generate only controller
    warning('Main interface generated through controller generateVHDL method');
    ctrl = object.getController();
    circuit_parameters.generateInterface = true;
    ctrl.generateVHDL(circuit_parameters);
end

