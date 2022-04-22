function explicitCtrlGenerateC(object, circuit_parameters)
%EXPLICITCTRLGENERATEC
% Generates C code for exact or approximate explicit MPC controller
% placed within an embedded system (observer + controller)
% 
% This is a private method.

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
nd = object.getNumberOfUnmeasurableInputs;
ny = object.getNumberOfOutputs;
nxref = numel(object.getController.getTrackingVariable);

range = object.range;

folder = circuit_parameters.folder;
initial_state = circuit_parameters.initialState;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE CONTROLLER FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctrlPar = circuit_parameters;
ctrlPar.generate_main = false;

ctrl = object.getController();
ctrlReport = ctrl.generateC(ctrlPar);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE OBSERVER FILE IF NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obsPar = circuit_parameters;
obsPar.range = object.range;

obsPar.generate_main = false;

obs = object.getObserver();
obsReport = obs.generateC(obsPar);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE MAIN FILE WITH OBSERVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = fopen([folder,'main.c'],'w');

if ~isa(object.dynSys,'pwaSys')

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
    fprintf(f,'#define nu_ES %d\n',nu);
    fprintf(f,'#define nx_ES %d\n',nx);
    fprintf(f,'#define np_ES %d\n',np);
    fprintf(f,'#define nd_ES %d\n',nd);
    fprintf(f,'#define ny_ES %d\n',ny);
    fprintf(f,'#define nxref_ES %d\n\n',nxref);


    fprintf(f,'void scaleInput(int xADC[np_ES+nxref_ES+ny_ES], float xscale[np_ES+nxref_ES+ny_ES]);\n');
    fprintf(f,'void scaleOutput(int uADC[nOut_OBS], float uscale[nOut_OBS]);\n');
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
    fprintf(f,'float state_predict[nOut_OBS] = {');
    for i=1:nx+nd-1
        fprintf(f,'%.10e,',initial_state(i));
    end
    fprintf(f,'%.10e};\n',initial_state(end));
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'float state_update[nOut_OBS] = {');
        for i=1:nx+nd-1
            fprintf(f,'%.10e,',initial_state(i));
        end
        fprintf(f,'%.10e};\n',initial_state(end));
    end

    fprintf(f,'float state_old[nOut_OBS] = {');
    for i=1:nx+nd-1
        fprintf(f,'%f,',initial_state(i));
    end
    fprintf(f,'%f};\n',initial_state(end));
    fprintf(f,'int inputVectorADC[np_ES+nxref_ES+ny_ES];//input vector,  [u0,...,un,p1,...,pm,y1,...,yq] \n');
    fprintf(f,'float inputVector[np_ES+nxref_ES+ny_ES];\n');
    fprintf(f,'float signalToObs[nIn_OBS];\n');
    fprintf(f,'float signalToCtrl[nDim_CTRL];\n');
    fprintf(f,'float u[nY_CTRL];//vector containing control value \n');
    fprintf(f,'int uADC[nY_CTRL];//vector containing control value \n\n');
    fprintf(f,'configADC();//config ADC \n');
    fprintf(f,'configDAC();//config DAC  \n\n');
    fprintf(f,'configTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',obs.getSamplingTime);
    fprintf(f,'while(1)\n{\n');
    fprintf(f,'while(!timerTic);//wait for interrupt activation\n');
    fprintf(f,'timerTic = 0;\n');
    fprintf(f,'readInputADC(inputVectorADC);\n');
    fprintf(f,'scaleInput(inputVectorADC,inputVector);\n');
    fprintf(f,'/* y = [y1,y2,...,yn,p1,..,pq] */\n\n');

    % generate signal for observer
    fprintf(f,'//fill observer input vector\n');
    fprintf(f,'for(i=0;i<np_ES;i++)\n\tsignalToObs[i+nu_ES] = inputVector[i];\n');
    fprintf(f,'for(i=0;i<ny_ES;i++)\n\tsignalToObs[i+nu_ES+np_ES] = inputVector[i+np_ES+nxref_ES];\n\n');

    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'update(signalToObs,state_predict,state_update);//update x with the kalmann observer\n');
    end

    fprintf(f,'observation = observation+1;\n');
    fprintf(f,'if(observation = control_predict_ratio)\n{\n');
    % generate signal for controller
    fprintf(f,'//fill controller input vector\n');
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'for(i=0;i<nx_ES;i++)\n\tsignalToCtrl[i] = state_update[i];\n\n');
    else
        fprintf(f,'for(i=0;i<nx_ES;i++)\n\tsignalToCtrl[i] = state_predict[i];\n\n');
    end;
    fprintf(f,'for(i=0;i<np_ES;i++)\n\tsignalToCtrl[i+nx_ES] = inputVector[i];\n\n');
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'for(i=0;i<nd_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES] = state_update[i+nx_ES];\n\n');
    else
        fprintf(f,'for(i=0;i<nd_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES] = state_predict[i+nx_ES];\n\n');
    end;
    fprintf(f,'for(i=0;i<nxref_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES+nd_ES] = inputVector[i+nu_ES+ny_ES];\n\n');

    fprintf(f,'control(signalToCtrl,u);//compute control\n');
    fprintf(f,'observation = 0;\n');
    fprintf(f,'scaleOutput(uADC,u);\n');
    fprintf(f,'setOutputDAC(uADC);\n}\n');
    % generate signal for observer
    fprintf(f,'//fill observer input vector\n');
    fprintf(f,'for(i=0;i<nu_ES;i++)\n\tsignalToObs[i] = u[i];\n\n');
    fprintf(f,'predict(signalToObs,state_old,state_predict);//predict x with the kalmann observer;\n');
    fprintf(f,'for(i=0;i<nOut_OBS;i++)\n\tstate_old[i] = state_predict[i];\n\n}\n');
    fprintf(f,'}\n\n');

    fprintf(f,'void scaleInput(int xADC[np_ES+nxref_ES+ny_ES], float xscale[np_ES+nxref_ES+ny_ES])\n{\n');
    for i=1:np
        deltaReal =  range.pmax(i)-range.pmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.pmax(i)+range.pmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;

        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i-1,deltaReal/deltaADC,i-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:nxref
        deltaReal =  range.xref(i)-range.xref(i);
        deltaADC = circuit_parameters.inputRange.max(i+np)-circuit_parameters.inputRange.min(i+np);
        meanReal = (range.xref(i)+range.xref(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+np)+circuit_parameters.inputRange.min(i+np))/2;

        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i+np-1,deltaReal/deltaADC,i-1+np,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:ny
        deltaReal =  range.ymax(i)-range.ymin(i);
        deltaADC = circuit_parameters.inputRange.max(i+np+nxref)-circuit_parameters.inputRange.min(i+np+nxref);
        meanReal = (range.ymax(i)+range.ymin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+np+nxref)+circuit_parameters.inputRange.min(i+np+nxref))/2;

        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i+np+nxref-1,deltaReal/deltaADC,i-1+np+nxref,(meanReal-deltaReal/deltaADC*meanADC));
    end
    fprintf(f,'}\n');

    fprintf(f,'void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL])\n{\n');
    for i=1:nu
        deltaReal =  range.umax(i)-range.umin(i);
        deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
        meanReal = (range.umax(i)+range.umin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;

        fprintf(f,'uADC[%d] = %.10e*uscale[%d]+%.10e;\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    fclose(f);

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
    memory_size = ctrlReport.memory_size+obsReport.memory_size;
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])





else

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
    fprintf(f,'#define nu_ES %d\n',nu);
    fprintf(f,'#define nx_ES %d\n',nx);
    fprintf(f,'#define np_ES %d\n',np);
    fprintf(f,'#define nd_ES %d\n',nd);
    fprintf(f,'#define ny_ES %d\n',ny);
    fprintf(f,'#define nxref_ES %d\n\n',nxref);


    fprintf(f,'void scaleInput(int xADC[np_ES+nxref_ES+ny_ES], float xscale[np_ES+nxref_ES+ny_ES]);\n');
    fprintf(f,'void scaleOutput(int uADC[nOut_OBS], float uscale[nOut_OBS]);\n');
    fprintf(f,'int findDynamic(float inputVect[np_ES+ny_ES+nxref_ES], float state[nx_ES+nd_ES]);\n');
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
    fprintf(f,'float state_predict[nOut_OBS] = {');
    for i=1:nx+nd-1
        fprintf(f,'%f,',initial_state(i));
    end
    fprintf(f,'%f};\n',initial_state(end));
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'float state_update[nOut_OBS] = {');
        for i=1:nx+nd-1
            fprintf(f,'%f,',initial_state(i));
        end
        fprintf(f,'%f};\n',initial_state(end));
    end

    fprintf(f,'float state_old[nOut_OBS] = {');
    for i=1:nx+nd-1
        fprintf(f,'%f,',initial_state(i));
    end
    fprintf(f,'%f};\n',initial_state(end));
    fprintf(f,'int inputVectorADC[np_ES+nxref_ES+ny_ES];//input vector,  [u0,...,un,p1,...,pm,y1,...,yq] \n');
    fprintf(f,'float inputVector[np_ES+nxref_ES+ny_ES];\n');
    fprintf(f,'float signalToObs[nIn_OBS];\n');
    fprintf(f,'float signalToCtrl[nDim_CTRL];\n');
    fprintf(f,'float u[nY_CTRL];//vector containing control value \n');
    fprintf(f,'int uADC[nY_CTRL];//vector containing control value \n');
    fprintf(f,'int i;\n\n');
    fprintf(f,'configADC();//config ADC \n');
    fprintf(f,'configDAC();//config DAC  \n\n');
    fprintf(f,'configTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',obs.getSamplingTime);
    fprintf(f,'while(1)\n{\n');
    fprintf(f,'while(!timerTic);//wait for interrupt activation\n');
    fprintf(f,'timerTic = 0;\n');
    fprintf(f,'readInputADC(inputVectorADC);\n');
    fprintf(f,'scaleInput(inputVectorADC,inputVector);\n');
    fprintf(f,'/* y = [y1,y2,...,yn,p1,..,pq] */\n\n');

    % generate signal for observer
    fprintf(f,'//fill observer input vector\n');
    fprintf(f,'for(i=0;i<np_ES;i++)\n\tsignalToObs[i+nu_ES] = inputVector[i];\n');
    fprintf(f,'for(i=0;i<ny_ES;i++)\n\tsignalToObs[i+nu_ES+np_ES] = inputVector[i+np_ES+nxref_ES];\n\n');
    fprintf(f,'i = findDynamic(inputVector,state_predict);\n');
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'update(signalToObs,state_predict,state_update,i);//update x with the kalmann observer\n');
        fprintf(f,'i = findDynamic(inputVector,state_update);\n');
    end

    fprintf(f,'observation = observation+1;\n');
    fprintf(f,'if(observation = control_predict_ratio)\n{\n');
    % generate signal for controller
    fprintf(f,'//fill controller input vector\n');
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'for(i=0;i<nx_ES;i++)\n\tsignalToCtrl[i] = state_update[i];\n\n');
    else
        fprintf(f,'for(i=0;i<nx_ES;i++)\n\tsignalToCtrl[i] = state_predict[i];\n\n');
    end;
    fprintf(f,'for(i=0;i<np_ES;i++)\n\tsignalToCtrl[i+nx_ES] = inputVector[i];\n\n');
    if isa(obs,'kalmanFilter') || isa(obs,'switchedKalmanFilter')
        fprintf(f,'for(i=0;i<nd_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES] = state_update[i+nx_ES];\n\n');
    else
        fprintf(f,'for(i=0;i<nd_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES] = state_predict[i+nx_ES];\n\n');
    end;
    fprintf(f,'for(i=0;i<nxref_ES;i++)\n\tsignalToCtrl[i+nx_ES+np_ES+nd_ES] = inputVector[i+nu_ES+ny_ES];\n\n');

    fprintf(f,'control(signalToCtrl,u,i);//compute control\n');
    fprintf(f,'observation = 0;\n');
    fprintf(f,'scaleOutput(uADC,u);\n');
    fprintf(f,'setOutputDAC(uADC);\n}\n');
    % generate signal for observer
    fprintf(f,'//fill observer input vector\n');
    fprintf(f,'for(i=0;i<nu_ES;i++)\n\tsignalToObs[i] = u[i];\n\n');
    fprintf(f,'predict(signalToObs,state_old,state_predict,i);//predict x with the kalmann observer;\n');
    fprintf(f,'for(i=0;i<nOut_OBS;i++)\n\tstate_old[i] = state_predict[i];\n\n}\n');
    fprintf(f,'}\n\n');

    fprintf(f,'void scaleInput(int xADC[np_ES+nxref_ES+ny_ES], float xscale[np_ES+nxref_ES+ny_ES])\n{\n');
    for i=1:np
        deltaReal =  range.pmax(i)-range.pmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.pmax(i)+range.pmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;

        fprintf(f,'xscale[%d] = %f*xADC[%d]+%f;\n',i-1,deltaReal/deltaADC,i-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:nxref
        deltaReal =  range.xref(i)-range.xref(i);
        deltaADC = circuit_parameters.inputRange.max(i+np)-circuit_parameters.inputRange.min(i+np);
        meanReal = (range.xref(i)+range.xref(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+np)+circuit_parameters.inputRange.min(i+np))/2;

        fprintf(f,'xscale[%d] = %f*xADC[%d]+%f;\n',i+np-1,deltaReal/deltaADC,i-1+np,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:ny
        deltaReal =  range.ymax(i)-range.ymin(i);
        deltaADC = circuit_parameters.inputRange.max(i+np+nxref)-circuit_parameters.inputRange.min(i+np+nxref);
        meanReal = (range.ymax(i)+range.ymin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+np+nxref)+circuit_parameters.inputRange.min(i+np+nxref))/2;

        fprintf(f,'xscale[%d] = %f*xADC[%d]+%f;\n',i+np+nxref-1,deltaReal/deltaADC,i-1+np+nxref,(meanReal-deltaReal/deltaADC*meanADC));
    end
    fprintf(f,'}\n');

    fprintf(f,'void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL])\n{\n');
    for i=1:nu
        deltaReal =  range.umax(i)-range.umin(i);
        deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
        meanReal = (range.umax(i)+range.umin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;

        fprintf(f,'uADC[%d] = %f*uscale[%d]+%f;\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    % generate findDynamic
    ff = object.controller.getInformation.controllers;
    nDyn = numel(ff);
    Hd = [eye(object.nx+object.np+object.nd);-eye(object.nx+object.np+object.nd)];
    Kd = [range.xmax(:); range.pmax(:); range.dmax(:); - range.xmin(:); -range.pmin(:); -range.dmin(:)];

    HH = [];
    KK = [];

    reg = cell(nDyn,1);

    for i=1:nDyn
        reg{i} = Polyhedron([Hd;ff(i).H],[Kd;ff(i).K]);
        HH = [HH; ff(i).H];
        KK = [KK; ff(i).K];
    end



    MOBYDICpars = getMOBYDICpars;
    tol = MOBYDICpars.roundtol;


    HK = [HH KK];

    % remove same edge
    HK = unique(HK,'rows');

    % Create edge set
    edge_set = 1:size(HK,1);

    % Remove domain edges from edge set
    HKd = [Hd Kd];
    HK = tol*round(HK/tol);
    HKd = tol*round(HKd/tol);
    rem = ismember(HK,HKd,'rows');
    edge_set(rem) = [];



    % remove same edge in opposite description
    i = 1;
    stop = 0;
    while ~stop
        ii = edge_set(i);
        currentHK = HK(ii,:);
        kk = find(ismember(HK,-currentHK,'rows'));
        jj = find(ismember(edge_set,kk));
        edge_set(jj) = [];

        if i == numel(edge_set)
            stop = 1;
        end
        i = i+1;
    end

    HKred = HK(edge_set,:);
    Hred = HKred(:,1:end-1);
    Kred = HKred(:,end);

    res = zeros(numel(Kred),nDyn);

    for i=1:nDyn
        centr = reg{i}.chebyCenter;
        if centr.r > 0
            res(:,i) = (Hred*centr.x <= Kred);
        else
            error('Empty dynamic found! Chech H and K matrix');
        end
    end
    fprintf(f,'int findDynamic(float inputVect[np_ES+ny_ES+nxref_ES], float state[nx_ES+nd_ES])\n{\n');
    fprintf(f,'const float Hedge[%d][nx_ES+np_ES+nd_ES] = {{',numel(Kred));
    for i=1:numel(Kred)-1
        for j=1:size(Hred,2)-1
            fprintf(f,'%f,',Hred(i,j));
        end
        fprintf(f,'%f},\n{',Hred(i,end));
    end
    for j=1:size(Hred,2)-1
        fprintf(f,'%f,',Hred(end,j));
    end
    fprintf(f,'%f}};\n\n',Hred(end,end));

    fprintf(f,'const float Kedge[%d] = {',numel(Kred));
    for i=1:numel(Kred)-1
        fprintf(f,'%f,\n',Kred(i));
    end
    fprintf(f,'%f};\n\n',Kred(end));
    fprintf(f,'int i,j,res;\n');
    fprintf(f,'float tmp;\n');
    fprintf(f,'res = 0;\n');
    fprintf(f,'for(i=0;i<%d;i++)\n{\n',numel(Kred));
    fprintf(f,'tmp = 0;\n');
    fprintf(f,'for(j=0;j<nx_ES;j++)\n');
    fprintf(f,'    tmp += Hedge[i][j]*state[j];\n');
    fprintf(f,'for(j=0;j<np_ES;j++)\n');
    fprintf(f,'    tmp += Hedge[i][j+nx_ES]*inputVect[j];\n');
    fprintf(f,'for(j=0;j<nd_ES;j++)\n');
    fprintf(f,'    tmp += Hedge[i][j+nx_ES+np_ES]*state[j+nx_ES];\n');
    fprintf(f,'tmp -= Kedge[i];\n');
    fprintf(f,'if (tmp <= 0)\n');
    fprintf(f,'    res += (1<<i);\n');
    fprintf(f,'}\n');
    intRes = zeros(1,nDyn);

    for i=1:nDyn
        rr = res(:,i);
        tmp = 0;
        for j=1:numel(Kred)
            tmp = tmp + rr(j)*2^(j-1);
        end
        intRes(i) = tmp;
    end

    fprintf(f,'if(res == %d)\n',intRes(1));
    fprintf(f,'    return 0;\n');
    for i=2:numel(intRes)
        fprintf(f,'else if(res == %d)\n',intRes(i));
        fprintf(f,'    return %d;\n',i-1);
    end
    fprintf(f,'else\n');
    fprintf(f,'   return 0;\n');
    fprintf(f,'}\n');


    fclose(f);

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
    memory_size = ctrlReport.memory_size+obsReport.memory_size;
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])

end
end

