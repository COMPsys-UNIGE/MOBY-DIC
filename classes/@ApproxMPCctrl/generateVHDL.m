%% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the MPCctrl
%
% SYNTAX
%
% object = generateVHDL(object)
% object = generateVHDL(object,opts)
%
%  opts is a structure with the following fields:
%
% * inputResolution: number of bits used to code inputs (u,p,y) (default
%         inputResolution = 12)
% * inputRange: variation range of the input. It is a struct with field
%               min and max (default full range)
% * inputRepresentation: representation of the input (signed,unsigned)
%                       (default signed)
%
% * coeffResolution: number of bits used to code the coefficients in memory
%                   (default coeffResolution = 12)
%
% * outputResolution: number of bits used to code outputs (x,d) (default
%         outputResolution = 12)
% * outputRange: variation range of the output. It is a struct with field
%               min and max (default full range)
% * outputRepresentation: representation of the output (signed,unsigned)
%                       (default signed)
%
% * range : variable variation range; it is a struct with fields xmin,
%           xmax, umin, umax, pmin, pmax, dmin, dmax, ymin, ymax. It must
%           be provided
%
% * frequency: indicates the frequency (in Hz) at which the circuit will
%              work. Default frequency = 50000000 (50 MHz).
%
% * generateInterface : if true the function generate also the interface
% of the observer, the multiplyer bank and the timer to timing sample
%
% It is possible to specify further options for the synthesis process.
% options is a structure with the following fields:
%
% * folder: destination folder where the VHDL files are saved (default is
%           pwag_ser_circuit or pwag_par_circuit followed by a progressive
%           number).
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Alberto Oliveri (alberto.oliveri@unige.it)
%
% Copyright is with:
%
% * Copyright (C) 2011 University of Genoa, Italy.

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

function varargout = generateVHDL(object,varargin)

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

circuit_parameters = object.ApproxMPCctrlVHDLset(circuit_parameters);

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
    if strcmp(circuit_parameters.inputRepresentation,'signed');
        numi = decimal2signed(fix(dscale),nbitout,0);
    else
        numi = decimal2unsigned(fix(dscale),nbitout,0);
    end
    do_bin{i} = numi.bin;
end

% generate package
object.writePackage(circuit_parameters.outputRepresentation,nbit,nbit_coeff,nbitout,sampling_latency,do_bin,folder);
%generate main block
if strcmp(circuit_parameters.inputRepresentation,'signed')
    copyfile([getvhdlpath,'/ApproxMPCctrl/controllerCircuit_signed.vhd'],[folder,'/controllerCircuit.vhd']);
else
    copyfile([getvhdlpath,'/ApproxMPCctrl/controllerCircuit_unsigned.vhd'],[folder,'/controllerCircuit.vhd']);
end

for i=1:numel(idx)
    dscale = (do(i)-(range.umax(i)+range.umin(i))/2)*(outputRange.max(i)-outputRange.min(i))/(range.umax(i)-range.umin(i))+(outputRange.max(i)+outputRange.min(i))/2;
    if strcmp(circuit_parameters.inputRepresentation,'signed');
    numi = decimal2signed(fix(dscale),nbitout,0);
    else
        numi = decimal2unsigned(fix(dscale),nbitout,0);
    end
    do_bin{i} = numi.bin;
end

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
%     object.write
    copyfile([getvhdlpath,'/ApproxMPCctrl/mulBank.vhd'],[folder,'/mulBank.vhd']);
    
    % generate input shifter
%     object.writeInputShift(meanInADC,circuit_parameters);
    
    % generate output shifter
%     object.writeOutputShift(meanOutADC,circuit_parameters);
    
    
    
    
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
        fprintf(fout,'\t\t%s reference: [%f %f] --> x%d_ref: [%s %s]\n',object.xnames{ii(i)},range.xmin(i),range.xmax(i),ii(i),cirmin.bin,cirmax.bin);
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
