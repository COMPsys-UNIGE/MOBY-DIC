%% SYNTHESIZE
% Generates the VHDL files describing the digital circuit implementing
% the pwas virtual sensor
%
% SYNTAX
%
% object = synthesize(object,circuit_parameters,[options])
%
% The vs object is returned as output because it is modified inside the
% method. circuit_parameters is a structure with the following fields:
%
% * nbit: number of bits used to code inputs (default nbit = 12)
% * nbit_coeff: number of bits used to code the value of the function in 
%               the vertices of the dimplicial partition (default 
%               nbit_coeff  = nbit)
% * type: string defining the type of architecture to be employed,
%         parallel ('parallel') or serial ('serial'). Default type =
%         'serial'.
% * frequency: indicates the frequency (in Hz) at which the circuit will
%              work. Default frequency = 20000000 (20 MHz).
% * samplingTime: indicates the sampling time of the circuit, i.e. the time
%                 between one input acquisition and the following one. If
%                 not provided it is set to the minimum available.
% * z0: initial condition for the estimate z
%
% It is possible to specify further options for the synthesis process.
% options is a structure with the following fields:
%
% * folder: destination folder where the VHDL files are saved (default is
%           vs_ser_circuit or vs_par_circuit followed by a progressive
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
% * Copyright (C) 2012 University of Genoa, Italy.

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

function object = synthesize (object, circuit_parameters, options)

% Retrieve information from circuit_parameters structure
if ~isfield(circuit_parameters,'nbit')
    nbit = 12;
else
    nbit = circuit_parameters.nbit;
end

if ~isfield(circuit_parameters,'nbit_coeff')
    nbit_coeff = nbit;
else
    nbit_coeff = circuit_parameters.nbit_coeff;
end

if ~isfield(circuit_parameters,'type')
    type = 'serial';
else
    type = circuit_parameters.type;
end

if ~isfield(circuit_parameters,'frequency')
    frequency = 20000000;
else
    frequency = circuit_parameters.frequency;
end

if ~isfield(circuit_parameters,'z0')
    z0 = 0;
else
    z0 = circuit_parameters.z0;
end

if ~isfield(circuit_parameters,'samplingTime')
    samplingTime = [];
else
    samplingTime = circuit_parameters.samplingTime;
end

if nbit < 4
    error('You must use at least 4 bits');
end

if ~strcmpi(type,'serial') && ~strcmpi(type,'parallel')
    error('Unknown type of architecture! Allowed types: serial and parallel')
end

if ~exist('options','var')
    options = [];
end

if ~object.reducedComplexity 
    error('Only reduced complexity virtual sensors can be implemented in digital circuits');
end

% Number of partitions
np = object.getNumberOfPartitions();

if ~all(np == np(1))
    error('Only virtual sensors with the same number of partition per dimension can be synthesized');
end
np = np(1);

% Set or create folder in which to store the VHDL files

if isfield(options,'folder')
    folder = options.folder;
else
    folder = '';
end

if ~isempty(folder)
    if strcmp(folder(end),'\')
        folder(end)='/';
    elseif ~strcmp(folder(end),'/')
        folder=[folder,'/'];
    end
else
    if strcmpi(type,'serial')
        created = 0;
        i = 1;
        while ~created
            testfolder = [pwd,'/vs_ser_circuit_',num2str(i),'/'];
            if ~exist(testfolder,'dir')
                folder = testfolder;
                created = 1;
            else
                i = i+1;
            end
        end
    else
        created = 0;
        i = 1;
        while ~created
            testfolder = [pwd,'/vs_par_circuit_',num2str(i),'/'];
            if ~exist(testfolder,'dir')
                folder = testfolder;
                created = 1;
            else
                i = i+1;
            end
        end
    end
end

if ~exist(folder,'dir')
    mkdir(folder)
end

disp(['Destination folder: ', folder]);

% Synthesize pwas functions

fpwas = object.getFunction;

nfpwas = numel(fpwas);

current_dir = pwd;

% Indicates if data at current time instants are used
current = object.isCurrent();

% Retrieve minimum and maximum weight between all pwas functions
maxw = zeros(1,nfpwas);
minw = zeros(1,nfpwas);
for i = 1:nfpwas
    w = fpwas(i).getWeights();
    maxw(i) = max(w);
    minw(i) = min(w);
end
w_max = max(maxw);
w_min = min(minw);

% Detect domain for unmeasurable output z
if current
    D = fpwas(2).getDomain();
else
    D = fpwas(1).getDomain();
end

D = D(:,end);

% Boundaries for z
zmin = D(1);
zmax = D(2);

% Matrices used to bring the weights from [w_min, w_max] to [0, 2^nbit_coeff-1]
Aw = (2^nbit_coeff-1)/(w_max-w_min);
Bw = -Aw*w_min;

% Complete circuit_parameters structure
circuit_parameters.vs = 1;
circuit_parameters.Aw = Aw;
circuit_parameters.Bw = Bw;

% Loop on all pwas functions, to synthesize all of them
for ii = 1:nfpwas
    
    % Temporary folder for ii-th pwas function
    opts.folder = [folder,'/fpwas',num2str(ii)];
    
    % Synthesize ii-th pwas function
    fpwas(ii) = fpwas(ii).synthesize(circuit_parameters,opts);
    
    % Move to folder containing the VHDL files of ii-th function
    cd(opts.folder)
    
    % Rename file names and entity names by adding a progressive number to
    % the name itself; this allows to avoid name collisions in the circuit
    % containing all pwas functions
    files = dir;
    files(1:2) = [];

    nfiles = numel(files);

    filenames = cell(nfiles,1);

    k = 1;

    for i = 1:nfiles
        tmpname = files(i).name;
        if strfind(tmpname,'.vhd')
            filenames{k} = files(i).name(1:end-4);
            k = k+1;
        end
    end

    filenames = filenames(1:k-2);

    filenames = [filenames; 'Input'];

    nfiles = numel(filenames);
    
    
    for i = 1:nfiles
        
        fin = fopen([filenames{i},'.vhd'],'r');
        if fin ~= -1
            fout = fopen(['../',filenames{i},'_',num2str(ii),'.vhd'],'w');
            
            while ~feof(fin)
                s = fgetl(fin);
                for j = 1:nfiles
                    s = strrep(s, [filenames{j}], [filenames{j},'_',num2str(ii)]);
                end
                fprintf(fout,'%s\n',s);
            end
            
            fclose(fin);
            fclose(fout);
        end
        
    end
    
    cd(current_dir)
    
    rmdir(opts.folder,'s');

end

% Latency for input acquisition
acquisition_latency = 0;

% Number of inputs
nu = object.nu;
% Number of measured outputs
ny = object.ny;
% Time window for input
mu = object.mu;
% Time window for outputs
my = object.my;
% Autoregressive time window
mz = object.n;

% Build matrix T
% Each row of T corresponds to a time instant
% Each column corresponds to a variable, in this way:
% u_1 ... u_mu y_1 ... y_my z
% The matrix contains only zeros and ones, the opnes indicates the variable
% involved for each time instants.
% For example:
%
%     u_1 u_2 y_1 y_2 z
%  k   1   1   1   1  0
% k-1  1   1   1   1  1 
% k-2  1   1   0   0  0
%
% Means that z = f1(u_1(k), u_2(k), y_1(k), y_2(k)) + 
%              + f2(u_1(k-1), u_2(k-1), y_1(k-1), y_2(k-1), z(k-1)) + 
%              + f3(u_1(k-1), u_2(k-2)
T = zeros(numel(fpwas),nu+ny+1);
if current
    for i = 1:nu
        T(1:mu(i),i) = 1;
    end
    for i = 1:ny
        T(1:my(i),nu+i) = 1;
    end
    T(2:mz+1,end) = 1;
else
    for i = 1:nu
        T(1:mu(i),i) = 1;
    end
    for i = 1:ny
        T(1:my(i),nu+i) = 1;
    end
    T(1:mz,end) = 1;
end

% Construct matrices A and B for the scaling of inputs, by taking pieces of
% matrix from all pwas functions

tmpAu = cell(numel(fpwas),1);
tmpAy = cell(numel(fpwas),1);
tmpAz = cell(numel(fpwas),1);
tmpBu = cell(numel(fpwas),1);
tmpBy = cell(numel(fpwas),1);
tmpBz = cell(numel(fpwas),1);

latency = 0;

% Loop on all pwas functions
for i = 1:numel(fpwas)
    
    % synthesisInfo structure of the i-th pwas function
    tmpInfo = fpwas(i).getSynthesisInfo();
    
    % Take maximum latency
    if tmpInfo.latency > latency
        latency = tmpInfo.latency;
    end
    % A and B matrices for i-th pwas function
    tmpA = tmpInfo.A;
    tmpB = tmpInfo.B;
    
    % Number of inputs involved in the i-th pwas function
    numu = numel(find(T(i,1:nu)));
    
    % Extract elements of A and B involving u
    tmpAu{i} = tmpA(1:numu,1:numu);
    tmpBu{i} = tmpB(1:numu);
    % Remove extracted elements from original matrices
    tmpA(1:numu,:) = [];
    tmpA(:,1:numu) = [];
    tmpB(1:numu) = [];
    
    % Number of measured outputs involved in the i-th pwas function
    numy = numel(find(T(i,nu+1:nu+ny)));
    
    % Extract elements of A and B involving y
    tmpAy{i} = tmpA(1:numy,1:numy);
    tmpBy{i} = tmpB(1:numy);
    % Remove extracted elements from original matrices
    tmpA(1:numy,:) = [];
    tmpA(:,1:numy) = [];
    tmpB(1:numy) = [];
        
    % Number of unmeasured outputs involved in the i-th pwas function
    numz = numel(find(T(i,nu+ny+1:end)));
    
    % Extract elements of A and B involving y
    tmpAz{i} = tmpA(1:numz,1:numz);
    tmpBz{i} = tmpB(1:numz);
end

% Add 2 clock cycles to the latency
latency = latency+2;

% Sampling latency
if isempty(samplingTime)
    sampling_latency = latency+3;
else
    sampling_latency = round(samplingTime*frequency);
end

% Put together the pieces of matrices A and B to construct the global
% matrices
A = [];
B = [];
Az = [];
Bz = [];

for i = 1:numel(tmpAu)
    if ~isempty(tmpAu{i})
        A = blkdiag(A,tmpAu{i});
        B = [B;tmpBu{i}]; 
        break;
    end
end

for i = 1:numel(tmpAy)
    if ~isempty(tmpAy{i})
        A = blkdiag(A,tmpAy{i});
        B = [B;tmpBy{i}]; 
        break;
    end
end

for i = 1:numel(tmpAz)
    if ~isempty(tmpAz{i})
        Az = blkdiag(Az,tmpAz{i});
        Bz = [Bz;tmpBz{i}]; 
        break;
    end
end

% alpha and beta are the matrices used to bring the result z to the
% original range
alpha = 1/Az;
beta = -Bz/Az;

% Extract synthesisInfo structure from first pwas function
synthesisInfo = fpwas(1).getSynthesisInfo();

% Number of bits to code decimal part
ndec = synthesisInfo.nbit - synthesisInfo.nint;

% gamma and fmin are used into the circuit to bring the output z in the correct range
% in order to make the feedback (for dynamical virtual sensors)
% zin = (zout-fmin)/gamma
gamma = np*(2^ndec)/(Aw*(zmax-zmin));
fmin = Aw*zmin+nfpwas*Bw;

% Number of bits
nbit = synthesisInfo.nbit;
% Number of bits of integer part
nint = synthesisInfo.nint;
% Number ofbits of decimal part
ndec = nbit-nint;
% Number of bits 
nbit_coeff = synthesisInfo.nbit_coeff;

% Filling fields of structure synthesisInfo
object.synthesisInfo.folder = folder;
object.synthesisInfo.A = A;
object.synthesisInfo.B = B;
object.synthesisInfo.nbit = nbit;
object.synthesisInfo.nint = nint;
object.synthesisInfo.nbit_coeff = nbit_coeff;
object.synthesisInfo.nbitout = synthesisInfo.nbitout;
object.synthesisInfo.alpha = alpha;
object.synthesisInfo.beta = beta;
object.synthesisInfo.type = synthesisInfo.type;
object.synthesisInfo.frequency = synthesisInfo.frequency;
object.synthesisInfo.latency = synthesisInfo.latency;
object.synthesisInfo.sampling_latency = synthesisInfo.sampling_latency;
object.synthesisInfo.representation = 'unsigned';

% Write VHDL files
writePackage(object,acquisition_latency,sampling_latency,z0,gamma,fmin);
object.writeInterface(circuit_parameters,folder);
object.writeTest(synthesisInfo,folder);

% Number of bits for the output
nbitout = nbit;
% Number of bits for the integer part of the output
nintout = nint;
% Number of bits for the decimal part of the output
ndecout = ndec;

% Fill other fields of synthesisInfo
object.synthesisInfo.nbitout = synthesisInfo.nbit;
object.synthesisInfo.nintout = synthesisInfo.nint;

% Copy missing VHDL files
dirvhdl = getvhdlpath();
copyfile([dirvhdl,'/vs/FIFO_u.vhd'],[folder,'/FIFO_u.vhd']);
copyfile([dirvhdl,'/vs/FIFO_y.vhd'],[folder,'/FIFO_y.vhd']);
copyfile([dirvhdl,'/vs/FIFO_z.vhd'],[folder,'/FIFO_z.vhd']);

% Write log file
if strcmpi(type,'parallel')
    acquisition = 'parallel';
    memory_size = zeros(1,nfpwas);
    multipliers = zeros(1,nfpwas);
    for i = 1:nfpwas
        nploc = fpwas(i).getNumberOfPartitions();
        nfun = fpwas(i).getNumberOfFunctions();
        ndim = fpwas(i).getNumberOfDimensions();
        memory_size(1) = prod(nploc+1)*nfun*(ndim+1);
        multipliers = ndim+1;
    end
    architecture = 'vs_par';
else
    acquisition = 'parallel';
    memory_size = zeros(1,nfpwas);
    multipliers = zeros(1,nfpwas);
    for i = 1:nfpwas
        nploc = fpwas(i).getNumberOfPartitions();
        nfun = fpwas(i).getNumberOfFunctions();
        memory_size(i) = prod(nploc+1)*nfun;
        multipliers(i) = 1;
        architecture = 'vs_ser';
    end
end
memory_size = sum(memory_size);
multipliers = sum(multipliers);
    
object.synthesisInfo.multipliers = multipliers;
object.synthesisInfo.memsize = nbit_coeff*memory_size/8;

filename = strcat(folder,'synthesisInfo.log');

fout = fopen(filename, 'w');

fprintf(fout,'-------------------------------------------------------------\n');
fprintf(fout,'|                Circuit information summary                 |\n');
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Circuit architecture: %s\n',architecture);
fprintf(fout,'Input acquisition: %s\n\n',acquisition);
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Data representation:\n');
fprintf(fout,'\t - Representation = fixed point, unsigned\n');
fprintf(fout,'\t - Number of bits for the inputs:\n');
fprintf(fout,'\t\t - Total = %d\n',nbit);
fprintf(fout,'\t\t - Integer part = %d\n',nint);
fprintf(fout,'\t\t - Decimal part = %d\n',ndec);
fprintf(fout,'\t - Number of bits for the coefficients:\n');
fprintf(fout,'\t\t - Total = %d\n',nbit_coeff);
fprintf(fout,'\t\t - Integer part = %d\n',nbit_coeff);
fprintf(fout,'\t\t - Decimal part = %d\n',0);
fprintf(fout,'\t - Number of bits for the outputs:\n');
fprintf(fout,'\t\t - Total = %d\n',nbitout);
fprintf(fout,'\t\t - Integer part = %d\n',nintout);
fprintf(fout,'\t\t - Decimal part = %d\n\n',ndecout);
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Timings:\r\n');
if frequency < 1000
    fprintf(fout,'\t - Working frequency = %f Hz\n',frequency);
elseif frequency < 1000000
    fprintf(fout,'\t - Working frequency = %f kHz\n',frequency/1e3);
else
    fprintf(fout,'\t - Working frequency = %f MHz\n',frequency/1e6);
end

time = (latency+1)/frequency;
sampling_time = sampling_latency/frequency;

if time > 1
    fprintf(fout,'\t - Latency = %.2f s (%d clock cycles)\n',time,latency+1);
elseif time > 1e-3
    fprintf(fout,'\t - Latency = %.2f ms (%d clock cycles)\n',time*1e3,latency+1);
elseif time > 1e-6
    fprintf(fout,'\t - Latency = %.2f us (%d clock cycles)\n',time*1e6,latency+1);
else
    fprintf(fout,'\t - Latency = %.2f ns (%d clock cycles)\n',time*1e9,latency+1);
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

fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Resources:\n');
fprintf(fout,'\t - Multiplier(s) = %d\n\n',multipliers);
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'Memory size:\n');
fprintf(fout,'\t - Number of cells = %d\n',memory_size);
fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*memory_size/8);
fclose(fout);
disp(['See ', folder ,'synthesisInfo.log for implementation details'])

