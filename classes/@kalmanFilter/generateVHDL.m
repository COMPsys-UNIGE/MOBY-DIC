%% generateVHDL
% Generates the VHDL files describing the digital circuit implementing
% the kalman filter
%
% SYNTAX
%
% generateVHDL(object,circuit_parameters)
%
%  circuit_parameters is a structure with the following fields:
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
% * initialState: indicates the initial state of the observer. It must be a
%                 matrix build as follow:
%                   [x_init d_init]'
%                  Default value : [0 ... 0]'
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

function varargout = generateVHDL(object,circuit_parameters)

nx = object.getNumberOfStates;
nd = object.getNumberOfUnmeasurableInputs;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;

sampling_time = object.getSamplingTime;

circuit_parameters = observerVHDLset(object,circuit_parameters);

% Read circuit_parameters structure

nbit = circuit_parameters.inputResolution;

nbitout = circuit_parameters.outputResolution;

nbit_coeff = circuit_parameters.coeffResolution;

frequency = circuit_parameters.frequency;

initialState = circuit_parameters.initialState;

folder = circuit_parameters.folder;

range = circuit_parameters.range;

generateInterface = circuit_parameters.generateInterface;

inputRange = circuit_parameters.inputRange;

outputRange = circuit_parameters.outputRange;

deltaInADC = inputRange.max-inputRange.min;

meanInADC =  ceil((inputRange.max+inputRange.min)/2);

deltaOutADC = outputRange.max-outputRange.min;

meanOutADC =  ceil((outputRange.max+outputRange.min)/2);

deltaRealIn = [range.umax(:) - range.umin(:);...
    range.pmax(:) - range.pmin(:);...
    range.ymax(:) - range.ymin(:)]';

meanRealIn = [(range.umax(:)+range.umin(:))/2;...
    (range.pmax(:)+range.pmin(:))/2;...
    (range.ymax(:)+range.ymin(:))/2];

deltaRealOut = [range.xmax(:) - range.xmin(:);...
    range.dmax(:) - range.dmin(:)]';

meanRealOut = [(range.xmax(:)+range.xmin(:))/2;...
    (range.dmax(:)+range.dmin(:))/2];

mulScale =  deltaRealIn(:)./deltaInADC(:);
mulScale = mulScale(:)';

outScale = deltaRealOut(:)./deltaOutADC(:);
outScale = outScale(:)';

alphaZ = diag(outScale);
alphaQ = diag(mulScale);

matrScale = repmat(mulScale,nx+nd,1);
matrOutScale = repmat(outScale,nx+nd,1);

bound = [ (outputRange.min(:)-meanOutADC(:))' (inputRange.min(:)-meanInADC(:))';
    (outputRange.max(:)-meanOutADC(:))' (inputRange.max(:)-meanInADC(:))'];

sampling_latency= sampling_time * frequency;

disp(['Destination folder for VHDL files: ', folder]);

% get observer matrices
[Aobs, Bobs, Cobs, Dobs, Gxobs, Gyobs] = object.getMatrices();

Ascale = inv(alphaZ)*Aobs*alphaZ;
Bscale = inv(alphaZ)*Bobs*alphaQ;%.*repmat(mulScale,nx+nd,1)./repmat(outScale(:),1,np+ny+nu);
% Gxscale = Gxobs./outScale'+meanRealOut./outScale'-Bobs*meanRealIn./outScale'-Aobs*meanRealOut(:)./outScale';
Gxscale = -inv(alphaZ)*meanRealOut+inv(alphaZ)*Aobs*meanRealOut+inv(alphaZ)*Gxobs+inv(alphaZ)*Bobs*meanRealIn;

% Take only last rows of matrices
Cobs = Cobs(object.ny+1:end,:);
Dobs = Dobs(object.ny+1:end,:);
Gyobs = Gyobs(object.ny+1:end,:);

Cscale = inv(alphaZ)*Cobs*alphaZ;
Dscale = inv(alphaZ)*Dobs*alphaQ
Gyscale = -inv(alphaZ)*meanRealOut+inv(alphaZ)*Cobs*meanRealOut+inv(alphaZ)*Gyobs+inv(alphaZ)*Dobs*meanRealIn;

%Gyscale = Gyobs./outScale'+meanRealOut./outScale'-Dobs*meanRealIn./outScale'-Cobs*meanRealOut(:)./outScale';

% compute rescale matrix and compute bit for prediction matrices
[FintPred, GintPred, alphaPred, betaPred, ~, nint_coeffPred] = quantizeData([Ascale Bscale],Gxscale,bound,nbit,nbit_coeff);
% compute rescale matrix and compute bit for update matrices
[FintUpdate, GintUpdate, alphaUpdate, betaUpdate, nint, nint_coeffUpdate] = quantizeData([Cscale Dscale],Gyscale,bound,nbit,nbit_coeff);

disp('Generating VHDL files...');

if ~exist(folder,'dir')
    mkdir(folder)
end

initialState = inv(alphaZ)*(initialState-meanRealOut);

% generate package
object.generatePackage(folder,FintPred, GintPred, alphaPred, betaPred, nbit, nint, nbit_coeff, nint_coeffPred, FintUpdate, GintUpdate, alphaUpdate, betaUpdate, nint_coeffUpdate,initialState,sampling_latency);

% %generate main block
% object.writeFilter(folder);

% copy other file
copyfile([getvhdlpath,'/kalmanFilter/shifterMACPredictor.vhdl'],[folder,'/shifterMACPredictor.vhdl']);
copyfile([getvhdlpath,'/kalmanFilter/shifterMACUpdate.vhdl'],[folder,'/shifterMACUpdate.vhdl']);
copyfile([getvhdlpath,'/kalmanFilter/statePredict.vhdl'],[folder,'/statePredict.vhdl']);
copyfile([getvhdlpath,'/kalmanFilter/stateUpdate.vhdl'],[folder,'/stateUpdate.vhdl']);
copyfile([getvhdlpath,'/kalmanFilter/kalmanFilter.vhdl'],[folder,'/kalmanFilter.vhdl']);

memory_size = 2*(nx+nd)*(nx+nd+np+nu+ny+1);
    
    latency = nx+nd+np+nu+ny+2;
    Multiplier = (nx+nd);
    if nargout > 1
        error('Wrong number of outputs')
    elseif nargout == 1
        optO.latency = latency;
        optO.memory_size = memory_size;
        optO.multiplier = Multiplier;
        optO.range = range;
        optO.circuit_parameters = circuit_parameters;
        varargout{1} = optO;
    end

% if needed generate interface
if generateInterface
    object.writeFilterInterface(folder);
    copyfile([getvhdlpath,'/kalmanFilter/mulBank.vhdl'],[folder,'/mulBank.vhdl']);
    
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
    
    np = object.np;
    nu = object.nu;
    ny = object.ny;
    
    nd = object.nd;
    nx = object.nx;
    
    for i = 1:object.nu
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,cirmin.bin,cirmax.bin);
    end
    
    for i = 1:object.np
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nu),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nu),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+nu),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+nu),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
    end
    
    for i = 1:object.ny
        if strcmp(circuit_parameters.inputRepresentation,'signed')
            cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nu+np),nbit,0);
            cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nu+np),nbit,0);
        else
            cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+nu+np),nbit,0);
            cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i),nbit,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.ynames{i},range.ymin(i),range.ymax(i),i,cirmin.bin,cirmax.bin);
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Representation: %s\n',circuit_parameters.outputRepresentation);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    
    for i = 1:object.nx
        if strcmp(circuit_parameters.outputRepresentation,'signed')
            ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        else
            ucirmin = decimal2unsigned(circuit_parameters.outputRange.min(i),nbitout,0);
            ucirmax = decimal2unsigned(circuit_parameters.outputRange.max(i),nbitout,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,ucirmin.bin,ucirmax.bin);
    end
        for i = 1:object.nd
        if strcmp(circuit_parameters.outputRepresentation,'signed')
            ucirmin = decimal2signed(circuit_parameters.outputRange.min(i+nx),nbitout,0);
            ucirmax = decimal2signed(circuit_parameters.outputRange.max(i+nx),nbitout,0);
        else
            ucirmin = decimal2unsigned(circuit_parameters.outputRange.min(i+nx),nbitout,0);
            ucirmax = decimal2unsigned(circuit_parameters.outputRange.max(i+nx),nbitout,0);
        end
        fprintf(fout,'\t\t%s: [%f %f] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,ucirmin.bin,ucirmax.bin);
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
    
    fprintf(fout,'\nTIMINGS (for each operation):\r\n');
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
    fprintf(fout,'\t - Multiplier(s) = %d (%d x %d bits)\n\n',Multiplier,max(nbit,nbit_coeff),max(nbit,nbit_coeff));
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*memory_size/8);
    fclose(fout);
    edit([folder ,'VHDL_report.log'])
end




