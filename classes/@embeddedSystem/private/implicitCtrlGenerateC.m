function implicitCtrlGenerateC(object, circuit_parameters)
%IMPLICITCTRLGENERATEC
% Generates C code for implicit MPC controller placed within an embedded
% system (observer + controller)
% 
% This is a private method.

% Number of variables
nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
nd = object.getNumberOfUnmeasurableInputs;
ny = object.getNumberOfOutputs;
nxref = numel(object.getController.getTrackingVariable);

% Signals range
range = object.range;

% Folder of C files
folder = circuit_parameters.folder;

% Initial observer state
initial_state = circuit_parameters.initialState;

    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% GENERATE CONTROLLER FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctrlPar = circuit_parameters;
ctrlPar.generate_main = false;

ctrl = object.getController();
ctrl.generateC(ctrlPar);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%$$%% GENERATE OBSERVER FILES  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obsPar = circuit_parameters;
obsPar.range = object.range;

obsPar.generate_main = false;

obs = object.getObserver();
obs.generateC(obsPar);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GENERATE MAIN FILE (CONTROLLER + OBSERVER) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'main.c'],'w');

fprintf(f,'/* Auto generated code for implementation of embedded system.\n');
fprintf(f,'   Final user must implement methods configADC(), configDAC()\n');
fprintf(f,'   readInputADC(y), setOutputDAC(u) and configTimerInterrupt(To)\n');
fprintf(f,'   as described above:\n');
fprintf(f,'   configADC() - in this function ADC must be set properly\n');
fprintf(f,'   configDAC() - in this function DAC must be set properly\n');
fprintf(f,'   configTimerInterrupt(float interruptTic) - in this function \n');
fprintf(f,'         timer and input must be set properly in order to\n');
fprintf(f,'         activate interrupt every To seconds\n');
fprintf(f,'   readInputADC(float *y) - in this function the progrma must\n');
fprintf(f,'         read value of all system output\n');
fprintf(f,'   setOutputDAC() - in this function the progrma must write\n');
fprintf(f,'         to DAC controls value\n');
fprintf(f,'   It is also necessary to call function interruptTimerRoutine()\n');
fprintf(f,'   into the timer interrupt routine.                            */\n');

fprintf(f,'#include <stdio.h>\n');
fprintf(f,'#include "controller.h"\n');
fprintf(f,'#include "observer.h"\n\n');

fprintf(f,'#define nY %d\n',ny);

fprintf(f,'void scaleInput(int inADC[nP+nRef+nY], float inscale[nP+nRef+nY]);\n');
fprintf(f,'void scaleOutput(int uADC[nU], float uscale[nU]);\n');
fprintf(f,'\n\n');
fprintf(f,'char timerTic = 0;\n');
fprintf(f,'\n\n');
fprintf(f,'void interruptTimerRoutine()\n{\n');
fprintf(f,'timerTic = 1;\n');
fprintf(f,'}\n\n');
fprintf(f,'\n\n');
fprintf(f,'int main()\n{\n');
fprintf(f,'const int control_predict_ratio =%d;//nPredict=T_ctrl/T_obs\n',floor(ctrl.getSamplingTime/obs.getSamplingTime));
fprintf(f,'int observation = 0;//nPredict=T_ctrl/T_obs\n');

fprintf(f,'int i;\n');
fprintf(f,'float state_predict[nX+nD] = {');
for i=1:nx+nd-1
    fprintf(f,'%.10e,',initial_state(i));
end
fprintf(f,'%.10e};\n',initial_state(end));
if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
    fprintf(f,'float state_update[nX+nD] = {');
    for i=1:nx+nd-1
        fprintf(f,'%.10e,',initial_state(i));
    end
    fprintf(f,'%.10e};\n',initial_state(end));
end

fprintf(f,'float state_old[nX+nD] = {');
for i=1:nx+nd-1
    fprintf(f,'%f,',initial_state(i));
end
fprintf(f,'%f};\n',initial_state(end));

fprintf(f,'float state_old[nX+nD] = {');
for i=1:nx+nd-1
    fprintf(f,'%f,',initial_state(i));
end
fprintf(f,'%f};\n',initial_state(end));

fprintf(f,'int inADC[nP+nRef+nY];//input vector\n');
fprintf(f,'int inscale[nP+nRef+nY];//scaled input vector\n');

fprintf(f,'float signalToObs[nIn_OBS];\n');

fprintf(f,'float x[nX];//state\n');

if np ~= 0
    fprintf(f,'float p[nP];//parameter\n');
end
if nd ~= 0
    fprintf(f,'float d[nD];//unmeasurable input\n');
end
if nxref ~= 0
    fprintf(f,'float ref[nRef];//reference state\n');
end

fprintf(f,'float u[nU];//vector containing control value \n\n');
fprintf(f,'int uADC[nU]\n');

fprintf(f,'configADC();//config ADC \n');
fprintf(f,'configDAC();//config DAC  \n\n');
fprintf(f,'configTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',obs.getSamplingTime);
fprintf(f,'while(1)\n{\n');
fprintf(f,'while(!timerTic);//wait for interrupt activation\n');
fprintf(f,'timerTic = 0;\n');

fprintf(f,'readInputADC(inADC);\n');

fprintf(f,'scaleInput(inADC,inscale);\n');

% generate signal for observer
fprintf(f,'//fill observer input vector\n');
fprintf(f,'for(i=0;i<nP;i++)\n\tsignalToObs[i+nU] = inscale[i];\n');
fprintf(f,'for(i=0;i<nY;i++)\n\tsignalToObs[i+nU+nP] = inscale[i+nP+nRef];\n\n');

if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
    fprintf(f,'update(signalToObs,state_predict,state_update);//update x with the Kalman observer\n');
end

fprintf(f,'observation = observation+1;\n');
fprintf(f,'if(observation = control_predict_ratio)\n{\n');

% generate signal for controller
fprintf(f,'//fill controller input vectors\n');

if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
    fprintf(f,'for(i=0;i<nX;i++)\n\tx[i] = state_update[i];\n\n');
else
    fprintf(f,'for(i=0;i<nX;i++)\n\tx[i] = state_predict[i];\n\n');
end

if np ~= 0
    fprintf(f,'for(i=0;i<nP;i++)\n\tp[i] = inscale[i];\n\n');
end

if nd ~= 0
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'for(i=0;i<nD;i++)\n\td[i] = state_update[i+nX];\n\n');
    else
        fprintf(f,'for(i=0;i<nD;i++)\n\td[i] = state_predict[i+nX];\n\n');
    end
end

if nxref ~= 0
    fprintf(f,'for(i=0;i<nRef;i++)\n\tref[i] = inscale[i+nU+nY];\n\n');
end

if nxref ~= 0
    if np ~= 0 && nd ~= 0
        fprintf(f,'control(x,u,p,d,ref);\n');
    elseif np ~= 0
        fprintf(f,'control(x,u,p,ref);\n');
    elseif nd ~= 0
        fprintf(f,'control(x,u,d,ref);\n');
    else
        fprintf(f,'control(x,u,ref);\n');
    end
else
    if np ~= 0 && nd ~= 0
        fprintf(f,'control(x,u,p,d);\n');
    elseif np ~= 0
        fprintf(f,'control(x,u,p);\n');
    elseif nd ~= 0
        fprintf(f,'control(x,u,d);\n');
    else
        fprintf(f,'control(x,u);\n');
    end
end

fprintf(f,'observation = 0;\n');
fprintf(f,'scaleOutput(uADC,u);\n');
fprintf(f,'setOutputDAC(uADC);\n}\n');

% generate signal for observer
fprintf(f,'//fill observer input vector\n');
fprintf(f,'for(i=0;i<nU;i++)\n\tsignalToObs[i] = u[i];\n\n');
fprintf(f,'predict(signalToObs,state_old,state_predict);//predict x with the Kalman observer;\n');
fprintf(f,'for(i=0;i<nX+nD;i++)\n\tstate_old[i] = state_predict[i];\n\n}\n');
fprintf(f,'}\n\n');

fprintf(f,'void scaleInput(int inADC[nP+nRef+nY], float inscale[nP+nRef+nY])\n{\n');
for i=1:np
    deltaReal =  range.pmax(i)-range.pmin(i);
    deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
    meanReal = (range.pmax(i)+range.pmin(i))/2;
    meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;

    fprintf(f,'inscale[%d] = %.10e*inADC[%d]+%.10e;\n',i-1,deltaReal/deltaADC,i-1,(meanReal-deltaReal/deltaADC*meanADC));
end
for i=1:nxref
    deltaReal =  range.xrefmax(i)-range.xrefmin(i);
    deltaADC = circuit_parameters.inputRange.max(i+np)-circuit_parameters.inputRange.min(i+np);
    meanReal = (range.xrefmax(i)+range.xrefmin(i))/2;
    meanADC = (circuit_parameters.inputRange.max(i+np)+circuit_parameters.inputRange.min(i+np))/2;

    fprintf(f,'inscale[%d] = %.10e*inADC[%d]+%.10e;\n',i+np-1,deltaReal/deltaADC,i-1+np,(meanReal-deltaReal/deltaADC*meanADC));
end
for i=1:ny
    deltaReal =  range.ymax(i)-range.ymin(i);
    deltaADC = circuit_parameters.inputRange.max(i+np+nxref)-circuit_parameters.inputRange.min(i+np+nxref);
    meanReal = (range.ymax(i)+range.ymin(i))/2;
    meanADC = (circuit_parameters.inputRange.max(i+np+nxref)+circuit_parameters.inputRange.min(i+np+nxref))/2;

    fprintf(f,'inscale[%d] = %.10e*inADC[%d]+%.10e;\n',i+np+nxref-1,deltaReal/deltaADC,i-1+np+nxref,(meanReal-deltaReal/deltaADC*meanADC));
end
fprintf(f,'}\n');

fprintf(f,'void scaleOutput(int uADC[nU], float uscale[nU])\n{\n');
for i=1:nu
    deltaReal =  range.umax(i)-range.umin(i);
    deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
    meanReal = (range.umax(i)+range.umin(i))/2;
    meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;

    fprintf(f,'uADC[%d] = %.10e*uscale[%d]+%.10e;\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
end
fprintf(f,'}\n');
fclose(f);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = strcat(folder,'C_report.log');

fout = fopen(filename, 'w');

nbit = circuit_parameters.inputResolution;
nbitout = circuit_parameters.outputResolution;

fprintf(fout,'-------------------------------------------------------------\n');
fprintf(fout,'|                Circuit information report                  |\n');
fprintf(fout,'-------------------------------------------------------------\n\n');
fprintf(fout,'INPUTS\n');
fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.inputResolution);
fprintf(fout,'\t - Range (model --> circuit):\n');
for i = 1:np
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
    fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
end
for i = 1:nxref
    track = object.getController.getTrackingVariable;
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i+np),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i+np),nbit,0);
    fprintf(fout,'\t\t%s_ref: [%f %f] --> xref%d: [%s %s]\n',object.xnames{track(i)},range.xrefmin(i),range.xrefmax(i),i,cirmin.bin,cirmax.bin);
end
for i = 1:ny
    cirmin = decimal2signed(circuit_parameters.inputRange.min(i+np+nxref),nbit,0);
    cirmax = decimal2signed(circuit_parameters.inputRange.max(i+np+nxref),nbit,0);
    fprintf(fout,'\t\t%s: [%f %f] --> y%d: [%s %s]\n',object.ynames{i},range.ymin(i),range.ymax(i),i,cirmin.bin,cirmax.bin);
end
fprintf(fout,'\nOUTPUTS\n');
fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
fprintf(fout,'\t - Range (model --> circuit):\n');
for i = 1:nu
    ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
    ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
    fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
end

fclose(fout);
edit([folder ,'C_report.log'])

end

