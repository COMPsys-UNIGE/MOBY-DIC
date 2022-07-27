function varargout = generateVHDL(object,varargin)
% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the pwag function
%
% generateVHDL(OBJ)
% Generates the VHDL files for the circuit implementation on FPGA of all
% components of the possibly vector PWAG function with the default options.
% Several VHDL files are generated, the top-level block is pwagFunction.vhd
% A report is also generated showing the main circuit features (latency,
% multipliers, etc.). A fixed point data representation is used.
%
% generateVHDL(OBJ,IDX)
% Generates the VHDL files for the circuit implementation on FPGA of the
% components of the vector PWAG function indicated by IDX, with the default
% options. IDX is a vector containing the indices of the components to
% implement. A report is also generated showing the main circuit features
% (latency, multipliers, etc.).
%
% generateVHDL(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation on FPGA of all
% components of the possibly vector PWAG function with custom options.
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
%                        generated. If the pwagFunction is part of a 
%                        controller, the interface must not be generated. 
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
    idx = [];
    circuit_parameters = [];
elseif nargin == 2
    idx = [];
    circuit_parameters = varargin{1};
elseif nargin == 3
    idx = varargin{1};
    circuit_parameters = varargin{2};
else
    error('Wrong input arguments');
end

% Get domain and codomain dimensions
ndim = object.getDomainDimensions();
nfun = object.getCodomainDimensions();


% Set default value for idx
if isempty(idx)
    idx = 1:nfun;
end

circuit_parameters = object.pwagVHDLset(circuit_parameters,idx);

nx = object.getDomainDimensions;
nu = numel(idx);

[Hd Kd] = object.getDomain;
P = Polyhedron(Hd,Kd);

if ~P.isBounded
    error(['PWAG function domain must be bounded in order'...
        ' to generate circuit file']);
end
Dp = P.outerApprox();
DVert = Dp.V;

minV = min(DVert);
maxV = max(DVert);

range.xmin = minV;
range.xmax = maxV;

VV = object.getVertices;
V = [];
for i=1:numel(VV)
    V = [V; VV{i}];
end

evalV = object.eval(V');

umin = min(evalV');
umax = max(evalV');

range.umin = umin(idx);
range.umax = umax(idx);

umin = repmat(Inf,nfun,1);
umax = repmat(-Inf,nfun,1);

for i=1:object.getNumberOfRegions
    r = object.getRegions(i);
    Ri = Polyhedron(r.H,r.K);
    vertici = Ri.V;
    for j=1:nfun
        umin(j) = min([umin(j),r.F(j,:)*vertici'+r.G(j)]);
        umax(j) = max([umax(j),r.F(j,:)*vertici'+r.G(j)]);
    end
end

range.umin = umin(idx)';
range.umax = umax(idx)';

if isfield(circuit_parameters,'range')
    if isfield(circuit_parameters.range,'xmax')
        if numel(circuit_parameters.range.xmax) == 1
          range.xmax = repmat(circuit_parameters.range.xmax,ndim,1);
        else
          range.xmax = circuit_parameters.range.xmax;
        end ;
    end
    if isfield(circuit_parameters.range,'xmin')
        if numel(circuit_parameters.range.xmin) == 1
          range.xmin = repmat(circuit_parameters.range.xmin,ndim,1);
        else
          range.xmin = circuit_parameters.range.xmin;
        end ;
    end

    if isfield(circuit_parameters.range,'umax')
        if numel(circuit_parameters.range.umax) == 1
          range.umax = repmat(circuit_parameters.range.umax,nu,1);
        else
          range.umax = circuit_parameters.range.umax;
        end ;
    end
    if isfield(circuit_parameters.range,'umin')
        if numel(circuit_parameters.range.umin) == 1
          range.umin = repmat(circuit_parameters.range.umin,nu,1);
        else
          range.umin = circuit_parameters.range.umin;
        end ;
    end

end


% Read circuit_parameters structure

nbit = circuit_parameters.inputResolution;

nbit_coeff = circuit_parameters.coeffResolution;

frequency = circuit_parameters.frequency;

folder = circuit_parameters.folder;

generateInterface = circuit_parameters.generateInterface;

nbitout = circuit_parameters.outputResolution;

inputRange = circuit_parameters.inputRange;

outputRange = circuit_parameters.outputRange;

deltaInADC = inputRange.max-inputRange.min;
meanInADC =  ceil((inputRange.max+inputRange.min)/2);


deltaOutADC = outputRange.max-outputRange.min;
meanOutADC = ceil((outputRange.max+outputRange.min)/2);


deltaRealIn = range.xmax-range.xmin;
meanRealIn = (range.xmax(:)+range.xmin(:))/2;

deltaRealIn(deltaRealIn == 0) = 1;

deltaRealOut = range.umax-range.umin;
meanRealOut = (range.umax(:)+range.umin(:))/2;

deltaRealOut(deltaRealOut == 0) = 1;

inScale =  deltaRealIn(:)./deltaInADC(:);
inScale = inScale(:)';

alphaZ = diag(inScale);

outScale =  deltaRealOut(:)./deltaOutADC(:);
outScale = outScale(:)';

alphaQ = diag(outScale);

% sampling_latency= sampling_time * frequency;



% Rescale all regions the matrices


% for i=1:object.getNumberOfRegions
%     ri = object.getRegions(i);
%     regions(i).H = ri.H.*repmat(inMatrScale(1,:),size(ri.H,1),1);
%     regions(i).K = ri.K -ri.H*meanRealIn;
%     regions(i).G = (ri.G(1:nu,:)+ri.F(1:nu,:)*meanRealIn-meanRealOut)./outScale';
%     regions(i).F = (ri.F(1:nu,:).*repmat(inMatrScale(1,:),size(ri.F(1:nu,:),1),1))./outMatrScale;
% end
%
% scaledPwag = pwagFunction(regions);
%
%     scaledPwag = scaledPwag.computeTree;
%
%
% state = treeExplore(scaledPwag);

if isempty(object.getTree)
    object = object.computeTree;
end

state = treeExplore(object);

rem = [];
for i = 1:numel(state)
    if isempty(state(i).leaf)
        rem = [rem,i];
    end
end

state(rem) = [];

[H K] = object.getEdges;
[F G] = object.getFunctions;

[Hred,Kred,Fred,Gred,stateMem,FGi_size] = object.reduceForCircuit(H,K,F,G,state,idx);
nbitState = ceil(log2(numel(stateMem)+1));

bigMeanRealOutVect = [];
alphaQBigVect = [];
for i=1:numel(idx)
    alphaQBigVect = [alphaQBigVect, repmat(outScale(i),1,FGi_size(i))];
    bigMeanRealOutVect = [bigMeanRealOutVect;repmat(meanRealOut(i),FGi_size(i),1)];
end

alphaQBig = diag(alphaQBigVect);
invalphaQBig = inv(alphaQBig);

Hmem = Hred*alphaZ;
Kmem = Kred -Hred*meanRealIn;

Fmem =invalphaQBig*Fred*alphaZ;
Gmem = invalphaQBig*Gred+invalphaQBig*Fred*meanRealIn-invalphaQBig*bigMeanRealOutVect;
% Fmem = (Fred(1:nu,:).*repmat(inMatrScale(1,:),size(Fred(1:nu,:),1),1))./outMatrScale;
% Gmem = (Gred(1:nu,:)+Fred(1:nu,:)*meanRealIn-meanRealOut)./outScale';


% Rescale matrices

D = [inputRange.min(:)';
    inputRange.max(:)'];

% compute rescale matrix and compute bit
[Hint Kint Fint Gint alpha beta nintcoeff] = ...
    object.quantizeData(Hmem,Kmem,Fmem,Gmem,D,nbit,nbit_coeff);

nbitMem = ceil(log2(size([Hint Kint ; Fint Gint],1)));

disp(['Destination folder for VHDL files: ', folder]);
disp('Generating VHDL files...');

if ~exist(folder,'dir')
    mkdir(folder)
end

object.writeMemory(Hint,-Kint,Fint,Gint,nbitMem,nbit_coeff,folder);

object.writeFSM(stateMem,nbitMem,numel(idx),folder);

% generate package
object.writePackage(nbit,nbit_coeff,nbitout,nbitMem,nbitState,nintcoeff,alpha,beta,idx,folder);
%generate main block
copyfile([getvhdlpath,'/pwagFunction/pwagFunctionCircuit.vhd'],[folder,'/pwagFunctionCircuit.vhd']);
copyfile([getvhdlpath,'/pwagFunction/compute.vhd'],[folder,'/compute.vhd']);

memory_size = numel([Hmem Kmem;Fmem Gmem]);

latency = (object.getDomainDimensions+3)*(object.getTree.maxDepth+3+numel(idx));

Multiplier = object.getDomainDimensions;

if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    optO.latency = latency;
    optO.memory_size = memory_size;
    optO.multiplier = Multiplier;
    optO.range = range;
    varargout{1} = optO;
end


% if needed generate interface
if generateInterface
    object.writeInterface(circuit_parameters);
    copyfile([getvhdlpath,'/pwagFunction/mulBank.vhd'],[folder,'/mulBank.vhd']);
    
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
    for i = 1:ndim
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,cirmin.bin,cirmax.bin);
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Representation: %s\n',circuit_parameters.outputRepresentation);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    for i = idx
        if strcmp(circuit_parameters.outputRepresentation,'signed')
            ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        else
            ucirmin = decimal2unsigned(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2unsigned(circuit_parameters.outputRange.max(i),nbitout,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.ynames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
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
    time = (latency)/frequency;
    if time > 1
        fprintf(fout,'\t - Latency = %.2f s (%d clock cycles)\n',time,latency);
    elseif time > 1e-3
        fprintf(fout,'\t - Latency = %.2f ms (%d clock cycles)\n',time*1e3,latency);
    elseif time > 1e-6
        fprintf(fout,'\t - Latency = %.2f us (%d clock cycles)\n',time*1e6,latency);
    else
        fprintf(fout,'\t - Latency = %.2f ns (%d clock cycles)\n',time*1e9,latency);
    end
    fprintf(fout,'\n');
    fprintf(fout,'RESOURCES\n');
    fprintf(fout,'\t - Multiplier(s) = %d (%d x %d bits)\n\n',Multiplier,max(nbit,nbit_coeff),max(nbit,nbit_coeff));
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*memory_size/8);
    fclose(fout);
    edit([folder ,'VHDL_report.log'])
    
    
    
end
