function varargout = generateVHDL(object, varargin)
% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the virtual sensor
%
% generateVHDL(OBJ)
% Generates the VHDL files for the circuit implementation of the virtual 
% sensor with the default options. A report is also generated showing the 
% main circuit features (latency, multipliers, etc.).
%
% generateVHDL(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation of the virtual 
% sensor with custom options. OPTS is a  structure with the following fields:
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
%  -  frequency: clock frequency (in Hz) of the FPGA. Default value:
%                50000000 (50 MHz).
%  -  folder: folder where to create VHDL files. Default value: progressive
%             folder.
% A report is also generated showing the main circuit features (latency, 
% multipliers, etc.).
%
% PERF = generateVHDL(OBJ,...)
% If an output argument is provided, the report is not generated and a 
% structure PERF containing the circuit performances is returned. 
% PERF is a structure with the following fields:
%  - latency: circuit latency (in clock cycles)
%  - memory_size: number of elements in the circuit memory
%  - multiplier: number of multipliers 
%  - range: structure with fields xmin, xmax, umin, umax representing the
%           ranges of the function inputs and outputs. 
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Alberto Oliveri (alberto.oliveri@unige.it)
% * Tomaso Poggi (tpoggi@essbilbao.org)
%
% Copyright is with:
%
% * Copyright (C) 2010-2011 University of Genoa, Italy.

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

% Check input arguments
if nargin == 1
    opts = [];
elseif nargin == 2
    opts = varargin{1};
else
    error('Wrong input arguments');
end

% Check if the virual sensor has been identified
if ~object.isIdentified
    error('The virtual sensor has not been identified.');
end

% Extract pwasFunction object(s)
fun = object.getFunction();
% Number of functions
npwas = numel(fun);

% TO DO
% Se reduced complexity e mz>0, solo signed
opts = vsVHDLset(object,opts);

% Get data
folder = opts.folder;
mu = object.getInputTimeWindow();
my = object.getOutputTimeWindow();
mz = object.getAutoregressiveTimeWindow();
current = object.isCurrent();
domain = object.getDomain;

% Current folder
thisfolder = pwd;

% Retrieve maximum and minimum weight among all functions
w = [];
for i = 1:npwas
    w = [w; fun(i).getWeights()];
end
minw = min(w);
maxw = max(w);

% Set options for VHDL generation of pwas functions
opts.generateInterface = 1;

% Find marices A and b which transform the weights in the range specified
% by opts.outputRange
if object.reducedComplexity
    if mz > 0
        if opts.outputRange.max < max(maxw,domain.zmax) || ...
                opts.outputRange.min > min(minw,domain.zmin)
            error(['OPTS.outputRange must be greater than range [',...
                num2str(min(minw,domain.zmin)),' ',num2str(max(maxw,domain.zmax)),']']);
        end
        % Constants which bring the minimum and maximum weight to the
        % minimum and maximum range respectively. b must be equal to 0
        % since the estimation os given as the sum of pwas functions, and
        % if b coefficients were present they would be summed up 
        cmax = opts.outputRange.max/max(maxw,domain.zmax);
        cmin = opts.outputRange.min/min(minw,domain.zmin);
    else
        if opts.outputRange.max < maxw || opts.outputRange.min > minw
            error('Output range is too small with respect to weights');
        end
        cmax = opts.outputRange.max/maxw;
        cmin = opts.outputRange.min/minw;
    end
    
    % Do not consider coefficients lower than 1
    if abs(cmin) < 1
        cmin = inf;
    end
    if abs(cmax) < 1
        cmax = inf;
    end
    
    % Take the minimum constant in order to link the weights to the output
    % range
    A = min(cmax,cmin);
    b = 0;
    
else
    
    % If the sensor is full complexity, coefficient b can be different from
    % 0
    A = (opts.outputRange.max-opts.outputRange.min)/(maxw-minw);
    b = opts.outputRange.max-A*maxw;
    
end

% If the autoregressive window is > 0, compute the range for the past
% estimates in input to the virtual sensor. It must be coherent with the
% output range
if mz > 0
    rangemin = round(A*domain.zmin+b);
    rangemax = round(A*domain.zmax+b);   
    
    % Transform also the initial condition
    opts.initialCondition = A*opts.initialCondition+b;
end
    
% Save the original input and output ranges
inputRange = opts.inputRange;
outputRange = opts.outputRange;

% Loop on all pwas functions
for i = 1:npwas
   
    % Folder where the VHDL files of the i-th pwas function are created
    curfolder = [folder,'pwas',num2str(i)];
    opts.folder = curfolder;
    
    % If the autoregressive window is > 0, the input range of the past
    % estimates must be made coherent with the output range
    if mz > 0 
        % Replicate the ranges on all dimensions
        opts.inputRange.min = repmat(inputRange.min,fun(i).getDomainDimensions,1);
        opts.inputRange.max = repmat(inputRange.max,fun(i).getDomainDimensions,1);
        
        % Modify the range corresponding to the past estimation
        if current
            if i > 1 && i <= mz+1
                opts.inputRange.min(end-mz+1:end) = rangemin;
                opts.inputRange.max(end-mz+1:end) = rangemax;
            end
        else
            if i <= mz
                opts.inputRange.min(end-mz+1:end) = rangemin;
                opts.inputRange.max(end-mz+1:end) = rangemax;
            end
        end
    end  
    
    % Weights of the current function
    wcur = fun(i).getWeights();
    % Transform the current weights in the output range 
    opts.outputRange.min = round(A*min(wcur)+b);
    opts.outputRange.max = round(A*max(wcur)+b);
    
    % Generate the VHDL files of the current pwas function
    perf(i) = fun(i).generateVHDL(opts);
    
    % Move to folder containing the VHDL files of ii-th function
    cd(curfolder)
    
    % Rename file names and entity names by adding a progressive number to
    % the name itself; this allows to avoid name collisions in the circuit
    % containing all pwas functions
    files = dir;
    files(1:2) = [];
    nfiles = numel(files);
    filenames = cell(nfiles,1);

    k = 1;
    for j = 1:nfiles
        tmpname = files(j).name;
        if strfind(tmpname,'.vhd')
            filenames{k} = files(j).name(1:end-4);
            k = k+1;
        end
    end
    filenames = filenames(1:k-1);
    nfiles = numel(filenames);
        
    for j = 1:nfiles
        fin = fopen([filenames{j},'.vhd'],'r');
        if fin ~= -1
            fout = fopen(['../',filenames{j},'_',num2str(i),'.vhd'],'w');
            
            while ~feof(fin)
                s = fgetl(fin);
                for k = 1:nfiles
                    s = strrep(s, [filenames{k}], [filenames{k},'_',num2str(i)]);
                end
                fprintf(fout,'%s\n',s);
            end
            
            fclose(fin);
            fclose(fout);
        end
        
    end
    
    cd(thisfolder)
    
    try
        rmdir(curfolder,'s');
    catch
       warning('Unable to delete folder'); 
    end
end

% Write the interface and the package of the virtual sensor
opts.folder = folder;
object.writeInterface(opts);
object.writePackage(opts,perf);

% Path where the VHDL files are stored
dirvhdl = getvhdlpath();

% Copy buffers
if my ~= 0
    copyfile([dirvhdl,'/virtualSensor/bufferY.vhd'],...
        [folder,'/bufferY.vhd']);
end
if mu ~= 0
    copyfile([dirvhdl,'/virtualSensor/bufferU.vhd'],...
        [folder,'/bufferU.vhd']);
end
if mz ~= 0
    copyfile([dirvhdl,'/virtualSensor/bufferZ.vhd'],...
        [folder,'/bufferZ.vhd']);
end

latency = 0;
multipliers = 0;
memory_size = 0;
for i = 1:npwas
    multipliers = multipliers+perf(i).multipliers;
    latency = max(latency,perf(i).latency);
    memory_size = memory_size+perf(i).memory_size;
end
% TO DO
latency = latency + 0;

% TO DO
% Check su Ts e latency

nbit_mul = perf(1).nbitmul;

nbit = opts.inputResolution;
nbit_coeff = opts.coeffResolution;

% Prepare the output sructure
if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    range.umin = object.domain.umin;
    range.umax = object.domain.umax;
    range.ymin = object.domain.ymin;
    range.ymax = object.domain.ymax;
    if mz > 0
        range.zmin = min(minw,domain.zmin);
        range.zmax = max(maxw,domain.zmax);
    else
        range.zmin = minw;
        range.zmax = maxw;
    end
    cir_range.umin = inputRange.min(1);
    cir_range.umax = inputRange.max(1);
    cir_range.ymin = inputRange.min(1);
    cir_range.ymax = inputRange.max(1);
    if mz > 0
        cir_range.zmin = A*min(minw,domain.zmin)+b;
        cir_range.zmax = A*max(maxw,domain.zmax)+b;
    else
        cir_range.zmin = A*minw+b;
        cir_range.zmax = A*maxw+b;
    end
    optO.latency = latency;
    optO.memory_size = memory_size;
    optO.multipliers = multipliers;
    optO.nbitmul = nbit_mul;
    optO.range = range;
    optO.cir_range = cir_range;
    varargout{1} = optO;
end

filename = strcat(folder,'VHDL_report.log');

fout = fopen(filename, 'w');

fprintf(fout,'-------------------------------------------------------------\n');
fprintf(fout,'|                Circuit information report                  |\n');
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Circuit architecture: %s\n',opts.architecture);
fprintf(fout,'INPUTS\n');
fprintf(fout,'\t - Resolution: %d bits\n',opts.inputResolution);
fprintf(fout,'\t - Representation: %s\n',opts.inputRepresentation);
fprintf(fout,'\t - Range (model --> circuit):\n');
if strcmp(opts.inputRepresentation,'signed')
    cirmin = decimal2signed(inputRange.min(1),nbit,0);
    cirmax = decimal2signed(inputRange.max(1),nbit,0);
else
    cirmin = decimal2unsigned(inputRange.min(1),nbit,0);
    cirmax = decimal2unsigned(inputRange.max(1),nbit,0);
end
for i = 1:object.nu
    fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},...
        object.domain.umin(i),object.domain.umax(i),i,cirmin.bin,cirmax.bin);
end
for i = 1:object.ny
    fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.ynames{i},...
        object.domain.ymin(i),object.domain.ymax(i),i,cirmin.bin,cirmax.bin);
end

fprintf(fout,'OUTPUTS\n');
fprintf(fout,'\t - Resolution: %d bits\n',opts.outputResolution);
fprintf(fout,'\t - Representation: %s\n',opts.outputRepresentation);
fprintf(fout,'\t - Range (model --> circuit):\n');
if strcmp(opts.inputRepresentation,'signed')
    cirmin = decimal2signed(outputRange.min,nbit,0);
    cirmax = decimal2signed(outputRange.max,nbit,0);
else
    cirmin = decimal2unsigned(outputRange.min,nbit,0);
    cirmax = decimal2unsigned(outputRange.max,nbit,0);
end
if mz > 0
    fprintf(fout,'\t\t%s: [%f %f] --> z: [%s %s]\n',object.znames{1},...
        min(minw,domain.zmin),max(maxw,domain.zmax),cirmin.bin,cirmax.bin);
else
    fprintf(fout,'\t\t%s: [%f %f] --> z: [%s %s]\n',object.znames{1},...
        minw,maxw,cirmin.bin,cirmax.bin);
end
fprintf(fout,'COEFFICIENTS\n');
fprintf(fout,'\t - Resolution: %d bits\n',opts.coeffResolution);
fprintf(fout,'\n');
fprintf(fout,'Timings:\r\n');
if opts.frequency < 1000
    fprintf(fout,'\t - Workin1g frequency = %f Hz\n',opts.frequency);
elseif opts.frequency < 1000000
    fprintf(fout,'\t - Working frequency = %f kHz\n',opts.frequency/1e3);
else
    fprintf(fout,'\t - Working frequency = %f MHz\n',opts.frequency/1e6);
end

time = latency/opts.frequency;

if time > 1
    fprintf(fout,'\t - Latency = %.2f s (%d clock cycles)\n',time,latency);
elseif time > 1e-3
    fprintf(fout,'\t - Latency = %.2f ms (%d clock cycles)\n',time*1e3,latency);
elseif time > 1e-6
    fprintf(fout,'\t - Latency = %.2f us (%d clock cycles)\n',time*1e6,latency);
else
    fprintf(fout,'\t - Latency = %.2f ns (%d clock cycles)\n',time*1e9,latency);
end

Ts = object.Ts;

if Ts > 1
    fprintf(fout,'\t - Sampling time = %.2f s (%d clock cycles)\n',Ts,Ts*opts.frequency);
elseif Ts > 1e-3
    fprintf(fout,'\t - Sampling time = %.2f ms (%d clock cycles)\n',Ts*1e3,Ts*opts.frequency);
elseif Ts > 1e-6
    fprintf(fout,'\t - Sampling time = %.2f us (%d clock cycles)\n',Ts*1e6,Ts*opts.frequency);
else
    fprintf(fout,'\t - Sampling time = %.2f ns (%d clock cycles)\n',Ts*1e9,Ts*opts.frequency);
end

fprintf(fout,'\n');
fprintf(fout,'Resources:\n');
fprintf(fout,'\t - Multiplier(s) = %d (%d x %d bits)\n\n',multipliers,nbit_mul,nbit_mul);
fprintf(fout,'Memory size:\n');
fprintf(fout,'\t - Number of cells = %d\n',memory_size);
fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*memory_size/8);
fclose(fout);
edit([folder ,'VHDL_report.log'])
