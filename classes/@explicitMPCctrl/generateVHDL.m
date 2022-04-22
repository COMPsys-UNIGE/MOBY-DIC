function varargout = generateVHDL(object,varargin)
% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the MPC controller
%
% generateVHDL(OBJ)
% Generates the VHDL files for the circuit implementation on FPGA of the 
% explicit MPC controller with the default options.
% Several VHDL files are generated, the top-level block is controllerInterface.vhd
% A report is also generated showing the main circuit features (latency,
% multipliers, etc.). A fixed point data representation is used.
%
% generateVHDL(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation on FPGA of the 
% explicit MPC controller with custom options.
% OPTS is a structure with the following fields:
%  -  architecture: can be either 'fast' (default) or 'small'. The first
%                   choice produces a circuit with lower latency and bigger
%                   size, the second one viceversa.
%  -  inputResolution: number of bits used to represent the function
%                      inputs. In general, it corresponds to the resolution
%                      of the analog-to-digital converters (ADC). Default
%                      value: 12.
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
%  -  coeffResolution: number of bits used to represent the coefficients.
%                      This value is in general bounded by the resolution
%                      of the multipliers available on the FPGA. Default
%                      value: inputResolution.
%  -  outputResolution: number of bits used to represent the function
%                       outputs. In general, it corresponds to the resolution
%                       of the digital-to-analog converters (DAC). Default
%                       value: coeffResolution.
%  -  outputRepresentation: can be either 'signed' or 'unsigned'. This also
%                           depends on the DACs. Default value: 'signed'.
%  -  outputRange: structure with fields min and max, indicating the minimum
%                  and maximum value of the circuit outputs. Default: maximum
%                  representable range.
%  -  frequency: clock frequency (in Hz) of the FPGA. This information is 
%                used to calculate the circuit latency.  Default value:
%                50000000 (50 MHz).
%  -  generateInterface: flag indicating if a circuit interface must be
%                        generated. If the controller is part of an 
%                        embedded system, the interface must not be generated. 
%                        Default value: true.
%  -  folder: folder where VHDL files are created. Default value: progressive
%             folder.
% A report is also generated showing the main circuit features (latency,
% multipliers, etc.).
%
% PERF = generateVHDL(OBJ,...)
% If an output argument is provided, the report is not automatically opened 
% and a structure PERF containing some information about the circuit 
% implementation is returned. PERF is a structure with the following fields:
%  - latency: circuit latency (in clock cycles)
%  - memory_size: number of elements in the circuit memory
%  - multiplier: number of multipliers
%  - range: structure with the following fields:
%            - xmin, xmax: actual ranges of the function inputs
%            - xmin_cir, xmax_cir: ranges of the function inputs for the
%                                  circuit
%            - umin, umax: actual ranges of the function outputs
%            - umin_cir, umax_cir: ranges of the function outputs for the
%                                  circuit
%  - scale_x: structure with fields A and b, which are the matrices used to 
%             transform the function inputs from their actual range to the
%             circuit range: x_cir = A x + b
%  - scale_u: structure with fields A and b, which are the matrices used to 
%             transform the function outputs from their circuit range to
%             the actual range: u = A u_cir + b
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

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nxref = numel(object.getTrackingVariable);

% Get domain and codomain dimensions
nfun = object.getNumberOfInputs();

circuit_parameters = MPCctrlVHDLset(object,circuit_parameters);

idx = 1:nfun;

func = object.getFunction;

subopt = circuit_parameters;
subopt.generateInterface = false;

% Read circuit_parameters structure

nbit = circuit_parameters.inputResolution;

nbit_coeff = circuit_parameters.coeffResolution;

frequency = circuit_parameters.frequency;

folder = circuit_parameters.folder;

generateInterface = circuit_parameters.generateInterface;

nbitout = circuit_parameters.outputResolution;

inputRange = circuit_parameters.inputRange;

outputRange = circuit_parameters.outputRange;


meanInADC =  ceil((inputRange.max+inputRange.min)/2);

meanOutADC = ceil((outputRange.max+outputRange.min)/2);


sampling_time = object.getSamplingTime;

sampling_latency= sampling_time * frequency;





disp(['Destination folder for VHDL files: ', folder]);
disp('Generating VHDL files...');

if ~exist(folder,'dir')
    mkdir(folder)
end



% generate function and controller VHDL
funOut = func.generateVHDL(1:nfun,subopt);

range.xmin = funOut.range.xmin(1:nx);
range.xmax = funOut.range.xmax(1:nx);

range.pmin = funOut.range.xmin(nx+1:nx+np);
range.pmax = funOut.range.xmax(nx+1:nx+np);

range.dmin = funOut.range.xmin(nx+np+1:nx+np+nd);
range.dmax = funOut.range.xmax(nx+np+1:nx+np+nd);

range.xrefmin = funOut.range.xmin(nx+np+nd+1:end);
range.xrefmax = funOut.range.xmax(nx+np+nd+1:end);

range.umin = funOut.range.umin;
range.umax = funOut.range.umax;

do = circuit_parameters.defaultOutput;
do_bin = cell(numel(idx),1);

for i=1:numel(idx)
    dscale = (do(i)-(range.umax(i)+range.umin(i))/2)*(outputRange.max(i)-outputRange.min(i))/(range.umax(i)-range.umin(i))+(outputRange.max(i)+outputRange.min(i))/2;
%     if strcmp(circuit_parameters.inputRepresentation,'signed')
    if strcmp(circuit_parameters.outputRepresentation,'signed')
        numi = decimal2signed(fix(dscale),nbitout,0);
    else
        numi = decimal2unsigned(fix(dscale),nbitout,0);
    end
    do_bin{i} = numi.bin;
end

% generate package
object.writePackage(nbit,nbit_coeff,nbitout,sampling_latency,do_bin,folder);
%generate main block
copyfile([getvhdlpath,'/explicitMPCctrl/controllerCircuit.vhd'],[folder,'/controllerCircuit.vhd']);




if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    optO = funOut;
    optO.range = range;
    optO.circuit_parameters = circuit_parameters;
    varargout{1} = optO;
end

% if needed generate interface
if generateInterface
    object.writeInterface(circuit_parameters);
    copyfile([getvhdlpath,'/explicitMPCctrl/mulBank.vhd'],[folder,'/mulBank.vhd']);
    
    % generate input shifter
    object.writeInputShift(meanInADC,circuit_parameters);
    
    % generate output shifter
    object.writeOutputShift(meanOutADC,circuit_parameters);
    
    
    
    
    filename = strcat(folder,'VHDL_report.log');
    
    fout = fopen(filename, 'w');
    
    
    fprintf(fout,'-------------------------------------------------------------\n');
    fprintf(fout,'|                Circuit information report                  |\n');
    fprintf(fout,'-------------------------------------------------------------\n\n');
    fprintf(fout,'INPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.inputResolution);
    fprintf(fout,'\t - Representation: %s\n',circuit_parameters.inputRepresentation);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    
    nx = object.nx;
    np = object.np;
    nd = object.nd;
    nxref = numel(object.getTrackingVariable);
    nu = object.nu;
    
    for i = 1:object.nx
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,cirmin.bin,cirmax.bin);
    end
    
    for i = 1:np
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+nx),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+nx),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
    end
    
    for i = 1:nd
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+nx+np),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+nx+np),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,cirmin.bin,cirmax.bin);
    end
    
    ii = object.getTrackingVariable;
    for i = 1:nxref
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np+nd),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np+nd),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+nx+np+nd),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+nx+np+nd),nbit,0);
        end
        fprintf(fout,'\t\t%s reference: [%f %f] --> x%d_ref: [%s %s]\n',object.xnames{ii(i)},range.xmin(ii(i)),range.xmax(ii(i)),ii(i),cirmin.bin,cirmax.bin);
    end
    
    
    
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Representation: %s\n',circuit_parameters.outputRepresentation);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    for i = 1:nu
        if strcmp(circuit_parameters.outputRepresentation,'signed')
            ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        else
            ucirmin = decimal2unsigned(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2unsigned(circuit_parameters.outputRange.max(i),nbitout,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.unames{i},funOut.range.umin(i),funOut.range.umax(i),i,ucirmin.bin,ucirmax.bin);
    end
    fprintf(fout,'\nCOEFFICIENTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.coeffResolution);
    fprintf(fout,'\n');
    fprintf(fout,'FREQUENCY:\r\n');
    if frequency < 1000
        fprintf(fout,'\t - Working frequency = %f Hz\n',frequency);
    elseif frequency < 1000000
        fprintf(fout,'\t - Working frequency = %f kHz\n',frequency/1e3);
    else
        fprintf(fout,'\t - Working frequency = %f MHz\n',frequency/1e6);
    end
    
    fprintf(fout,'\nTIMINGS:\r\n');
    time = (funOut.latency)/frequency;
    
    
    if time > 1
        fprintf(fout,'\t - Latency = %.2f s (%d clock cycles)\n',time,funOut.latency);
    elseif time > 1e-3
        fprintf(fout,'\t - Latency = %.2f ms (%d clock cycles)\n',time*1e3,funOut.latency);
    elseif time > 1e-6
        fprintf(fout,'\t - Latency = %.2f us (%d clock cycles)\n',time*1e6,funOut.latency);
    else
        fprintf(fout,'\t - Latency = %.2f ns (%d clock cycles)\n',time*1e9,funOut.latency);
    end
    
    if sampling_time >= 1
        fprintf(fout,'\t - Sampling time = %.2f s (%d clock cycles)\n\n',sampling_time,sampling_latency);
    elseif sampling_time >= 1e-3
        fprintf(fout,'\t - Sampling time = %.2f ms (%d clock cycles)\n\n',sampling_time*1e3,sampling_latency);
    elseif sampling_time >= 1e-6
        fprintf(fout,'\t - Sampling time = %.2f us (%d clock cycles)\n\n',sampling_time*1e6,sampling_latency);
    else
        fprintf(fout,'\t - Sampling time = %.2f ns (%d clock cycles)\n\n',sampling_time*1e9,sampling_latency);
    end
    
    
    fprintf(fout,'\n');
    fprintf(fout,'RESOURCES\n');
    fprintf(fout,'\t - Multiplier(s) = %d (%d x %d bits)\n\n',funOut.multiplier,max(nbit,nbit_coeff),max(nbit,nbit_coeff));
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',funOut.memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*funOut.memory_size/8);
    fclose(fout);
    edit([folder ,'VHDL_report.log'])
    
    
    
end
