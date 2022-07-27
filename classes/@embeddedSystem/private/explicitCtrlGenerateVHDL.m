function out = explicitCtrlGenerateVHDL(object, circuit_parameters)
%EXPLICITCTRLGENERATEVHDL
% Generates VHDL code for exact or approximate explicit MPC controller
% placed within an embedded system (observer + controller)
% 
% This is a private method.

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
nd = object.getNumberOfUnmeasurableInputs;
ny = object.getNumberOfOutputs;
nxref = numel(object.getController.getTrackingVariable);

circuit_parameters = object.embeddedSystemVHDLset(circuit_parameters);

frequency = circuit_parameters.frequency;


meanIn = ceil((circuit_parameters.inputRange.max+circuit_parameters.inputRange.min)/2);
meanOut = ceil((circuit_parameters.outputRange.max+circuit_parameters.outputRange.min)/2);


Unbiased_in.max = circuit_parameters.inputRange.max-meanIn;
Unbiased_in.min = circuit_parameters.inputRange.min-meanIn;
Unbiased_out.max = circuit_parameters.outputRange.max-meanOut;
Unbiased_out.min = circuit_parameters.outputRange.min-meanOut;

ncoeff = circuit_parameters.coeffResolution;

% generate controller
% es input are [p xref y], for controller we need [x p d xref]
ctrlInputRangeMax = [repmat(2^(ncoeff-1)-1,1,nx)';...
Unbiased_in.max(1:np);...
repmat(2^(ncoeff-1)-1,1,nd)';...
Unbiased_in.max(np+1:np+nxref)];

ctrlInputRangeMin = [repmat(-(2^(ncoeff-1)),1,nx)';...
Unbiased_in.min(1:np);...
repmat(-(2^(ncoeff-1)),1,nd)';...
Unbiased_in.min(np+1:np+nxref)];



ctrlPar = circuit_parameters;
ctrlPar.inputResolution = ctrlPar.coeffResolution;
ctrlPar.inputRepresentation = 'signed';
ctrlPar.inputRange.max = ctrlInputRangeMax;
ctrlPar.inputRange.min = ctrlInputRangeMin;

ctrlPar.outputResolution = circuit_parameters.outputResolution;
ctrlPar.outputRepresentation = 'signed';
ctrlPar.outputRange = Unbiased_out;
ctrlPar.generateInterface = false;

ctrl = object.getController();
ctrlReport = ctrl.generateVHDL(ctrlPar);

% generate observer
% es input are [p xref y], for observer we need [u p y]
obsInputRangeMax = [Unbiased_out.max;...
Unbiased_in.max(1:np);...
Unbiased_in.max(np+nxref+1:np+nxref+ny)];

obsInputRangeMin = [Unbiased_out.min;...
Unbiased_in.min(1:np);...
Unbiased_in.min(np+nxref+1:np+nxref+ny)];

obsPar = circuit_parameters;
obsPar.inputResolution = circuit_parameters.coeffResolution;
obsPar.inputRepresentation = 'signed';
obsPar.inputRange.max = obsInputRangeMax;
obsPar.inputRange.min = obsInputRangeMin;

obsPar.outputResolution = ncoeff;
obsPar.outputRepresentation = 'signed';
obsPar.outputRange.max = repmat(2^(ncoeff-1)-1,1,nx+nd);
obsPar.outputRange.min = repmat(-(2^(ncoeff-1)),1,nx+nd);
obsPar.range = object.range;

obsPar.generateInterface = false;

obs = object.getObserver();
obsReport = obs.generateVHDL(obsPar);


nbit = circuit_parameters.inputResolution;
nbit_coeff = circuit_parameters.coeffResolution;
nbitout = circuit_parameters.outputResolution;

nMul_obs = obsReport.multiplier;
nMul_ctrl = ctrlReport.multiplier;

nMul = max(nMul_obs,nMul_ctrl);

obs_latency = object.getObserver.getSamplingTime*circuit_parameters.frequency;
ctrl_latency = object.getController.getSamplingTime*circuit_parameters.frequency;
folder = circuit_parameters.folder;

object.writeInputShift(meanIn,circuit_parameters);
object.writeOutputShift(meanOut,circuit_parameters);

nDyn = [];
nedge = 0;
if isa(object.dynSys,'pwaSys')
nDyn = object.dynSys.getNumberOfDynamics;
nedge = object.writeFindDynamic(circuit_parameters);
end

eff_bit = nbit;
if isa(ctrl,'ApproxMPCctrl')||isa(ctrl,'switchedApproxMPCctrl')
eff_bit = eff_bit+1;
end

nbitMul = max(nbit_coeff,eff_bit); 

object.writePackage(nbit,nbit_coeff,nbitout,nMul,nbitMul,nMul_ctrl,nMul_obs,obs_latency,ctrl_latency,nDyn,folder);
% generate mulbank
copyfile([getvhdlpath,'/embeddedSystem/mulBank.vhd'],[folder,'/mulBank.vhd']);
% generate timingFSM
if isa(object.dynSys,'pwaSys')
if isa(object.getObserver(),'switchedKalmanFilter')
    copyfile([getvhdlpath,'/embeddedSystem/pwaSys_timingFSM_filter.vhd'],[folder,'/timingFSM.vhd']);
else
    copyfile([getvhdlpath,'/embeddedSystem/pwaSys_timingFSM_predictor.vhd'],[folder,'/timingFSM.vhd']);
end
else
if isa(object.getObserver(),'kalmanFilter')
    copyfile([getvhdlpath,'/embeddedSystem/ltiSys_timingFSM_filter.vhd'],[folder,'/timingFSM.vhd']);
else
    copyfile([getvhdlpath,'/embeddedSystem/ltiSys_timingFSM_predictor.vhd'],[folder,'/timingFSM.vhd']);
end
end

range = object.range;

object.writeEmbeddedSystem(circuit_parameters);

% Write report
filename = strcat(folder,'VHDL_report.log');

fout = fopen(filename, 'w');

mul = max(ctrlReport.multiplier,obsReport.multiplier);
% if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
% latency = ctrlReport.latency+2*obsReport.latency+nedge+3;
% else
    latency = ctrlReport.latency+2*obsReport.latency+nedge+3;
% end
memory_size = ctrlReport.memory_size+obsReport.memory_size+nedge*(nx+nd+np+1);

out.latency = latency;
out.memory_size = memory_size;
out.multiplier = mul;
out.range = object.range;
out.circuit_parameters = circuit_parameters;


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
nxref = numel(object.getController.getTrackingVariable);
nu = object.nu;

for i = 1:np
if strcmp(circuit_parameters.inputRepresentation,'signed')
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
else
    cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i),nbit,0);
    cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i),nbit,0);
end
fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
end

for i = 1:object.ny
if strcmp(circuit_parameters.inputRepresentation,'signed')
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i+np),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i+np),nbit,0);
else
    cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+np),nbit,0);
    cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+np),nbit,0);
end
fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.ynames{i},range.ymin(i),range.ymax(i),i,cirmin.bin,cirmax.bin);
end

ii = object.getController.getTrackingVariable;
for i = 1:nxref
if strcmp(circuit_parameters.inputRepresentation,'signed')
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i+np+ny),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i+np+nd),nbit,0);
else
    cirmin = decimal2unsigned(circuit_parameters.inputRange.min(i+np+ny),nbit,0);
    cirmax = decimal2unsigned(circuit_parameters.inputRange.max(i+np+ny),nbit,0);
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
fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
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
sampling_latency = obs_latency;
sampling_time = obs_latency/frequency;
if sampling_time >= 1
fprintf(fout,'\t - Observer sampling time = %.2f s (%d clock cycles)\n',sampling_time,sampling_latency);
elseif sampling_time >= 1e-3
fprintf(fout,'\t - Observer sampling time = %.2f ms (%d clock cycles)\n',sampling_time*1e3,sampling_latency);
elseif sampling_time >= 1e-6
fprintf(fout,'\t - Observer sampling time = %.2f us (%d clock cycles)\n',sampling_time*1e6,sampling_latency);
else
fprintf(fout,'\t - Observer sampling time = %.2f ns (%d clock cycles)\n',sampling_time*1e9,sampling_latency);
end
sampling_latency = ctrl_latency;
sampling_time = ctrl_latency/frequency;
if sampling_time >= 1
fprintf(fout,'\t - Controller sampling time = %.2f s (%d clock cycles)\n\n',sampling_time,sampling_latency);
elseif sampling_time >= 1e-3
fprintf(fout,'\t - Controller sampling time = %.2f ms (%d clock cycles)\n\n',sampling_time*1e3,sampling_latency);
elseif sampling_time >= 1e-6
fprintf(fout,'\t - Controller sampling time = %.2f us (%d clock cycles)\n\n',sampling_time*1e6,sampling_latency);
else
fprintf(fout,'\t - Controller sampling time = %.2f ns (%d clock cycles)\n\n',sampling_time*1e9,sampling_latency);
end

fprintf(fout,'\n');
fprintf(fout,'RESOURCES\n');
fprintf(fout,'\t - Multiplier(s) = %d (%d x %d bits)\n\n',mul,nbitMul,max(nbit,nbit_coeff));

fprintf(fout,'MEMORY SIZE\n');
fprintf(fout,'\t - Number of cells = %d\n',memory_size);
fprintf(fout,'\t - Word size = %d bits\n',nbit_coeff);
fprintf(fout,'\t - Total occupation = %.3f bytes\n',nbit_coeff*memory_size/8);
fclose(fout);
edit([folder ,'VHDL_report.log'])
end

