%% generateC
% Generates the C files describing the digital circuit implementing
% the ApproxMPCctrl
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
% * generate_main : if true main file will be generated (default true)
% * folder: destination folder where the C files are saved
% * outputResolution: number of bits used to code output (u,p,y) (default
%         inputResolution = 12)
% * outputRange: variation range of the output. It is a struct with field
%               min and max (default full range)
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

function varargout = generateC(object,varargin)

if nargin == 1
    circuit_parameters = [];
elseif nargin == 2
    circuit_parameters = varargin{1};
else
    error('Wrong input arguments');
end

circuit_parameters = object.observerCset(circuit_parameters);

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;

nu = object.getNumberOfInputs;
ny = object.getNumberOfInputs;
range = circuit_parameters.range;

nIn = np+nu+ny;
nOut = nx+nd;

initial_state = circuit_parameters.initialState;

folder = circuit_parameters.folder;

if ~exist(folder)
    mkdir(folder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE MAIN FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if circuit_parameters.generate_main
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
    fprintf(f,'#include "observer.h"\n');
    fprintf(f,'void scaleInput(int xADC[nIn_OBS], float xscale[nIn_OBS]);\n');
    fprintf(f,'void scaleOutput(int uADC[nOut_OBS], float uscale[nOut_OBS]);\n');
    fprintf(f,'int findDynamic(float inputVect[nIn_OBS], float state[nOut_OBS]);\n');
    fprintf(f,'\n\n');
    fprintf(f,'char observeGo = 0;\n');
    fprintf(f,'\n\n');
    fprintf(f,'void interruptTimerRoutine()\n{\n');
    fprintf(f,'observeGo = 1;\n');
    fprintf(f,'}\n\n');
    fprintf(f,'\n\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'int inADC[nIn_OBS];\n');
    fprintf(f,'float inVector[nIn_OBS];//input vector\n');
    fprintf(f,'float state_predict[nOut_OBS] = {');
    for i=1:nx+nd-1
        fprintf(f,'%.10e,',initial_state(i));
    end
    fprintf(f,'%.10e};\n',initial_state(end));
    fprintf(f,'float old_state[nOut_OBS];\n');
    fprintf(f,'int outADC[nOut_OBS];\n');
    fprintf(f,'int i;\n');
    fprintf(f,'configADC();//config ADC \n');
    fprintf(f,'configDAC();//config DAC  \n\n');
    fprintf(f,'configTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',object.getSamplingTime);
    fprintf(f,'while(1)\n{\n');
    fprintf(f,'while(!observeGo);//wait for interrupt activation\n');
    fprintf(f,'observeGo = 0;\n');
    fprintf(f,'/* Read form ADC in inADC = [u0,...,un,p1,...,pm,y1,...,yq];*/\n');
    fprintf(f,'readInputADC(inADC);\n');
    fprintf(f,'scaleInput(inADC,inVector);\n');
    fprintf(f,'i = findDynamic(inVector,state_predict);\n');
    fprintf(f,'predict(inVector,old_state,state_predict,i);\n');
    fprintf(f,'scaleOutput(outADC,state_predict);\n');
    fprintf(f,'setOutputDAC(outADC);\n');
    fprintf(f,'for(i=0;i<nOut_OBS;i++)\n\told_state[i] = state_predict[i];\n}\n');
    fprintf(f,'}\n\n');
    
    fprintf(f,'void scaleInput(int xADC[nIn_OBS], float xscale[nIn_OBS])\n{\n');
    for i=1:nu
        deltaReal =  range.umax(i)-range.umin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.umax(i)+range.umin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        
        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i-1,deltaReal/deltaADC,i-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:np
        deltaReal =  range.pmax(i)-range.pmin(i);
        deltaADC = circuit_parameters.inputRange.max(i+nu)-circuit_parameters.inputRange.min(i+nu);
        meanReal = (range.pmax(i)+range.pmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+nu)+circuit_parameters.inputRange.min(i+nu))/2;
        
        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i+nu-1,deltaReal/deltaADC,i-1+nu+np,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:ny
        deltaReal =  range.ymax(i)-range.ymin(i);
        deltaADC = circuit_parameters.inputRange.max(i+nu+np)-circuit_parameters.inputRange.min(i+nu+np);
        meanReal = (range.ymax(i)+range.ymin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i+nu+np)+circuit_parameters.inputRange.min(i+nu+np))/2;
        
        fprintf(f,'xscale[%d] = %.10e*xADC[%d]+%.10e;\n',i+nu+np-1,deltaReal/deltaADC,i-1+nu+np,(meanReal-deltaReal/deltaADC*meanADC));
    end
    fprintf(f,'}\n');
    
    fprintf(f,'void scaleOutput(int uADC[nOut_OBS], float uscale[nOut_OBS])\n{\n');
    for i=1:nx
        deltaReal =  range.xmax(i)-range.xmin(i);
        deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
        meanReal = (range.xmax(i)+range.xmin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;
        
        fprintf(f,'uADC[%d] = %.10e*uscale[%d]+%.10e;\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    
    for i=1:nd
        deltaReal =  range.dmax(i)-range.dmin(i);
        deltaADC = circuit_parameters.outputRange.max(i+nx)-circuit_parameters.outputRange.min(i+nx);
        meanReal = (range.dmax(i)+range.dmin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i+nx)+circuit_parameters.outputRange.min(i+nx))/2;
        
        fprintf(f,'uADC[%d] = %.10e*uscale[%d]+%.10e;\n',i+nx-1,deltaADC/deltaReal,i+nx-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    
    % generate findDynamic
    ff = object.filters;
    
    Hd = [eye(object.nx+object.np+object.nd);-eye(object.nx+object.np+object.nd)];
    Kd = [range.xmax(:); range.pmax(:); range.dmax(:); - range.xmin(:); -range.pmin(:); -range.dmin(:)];
    
    HH = [];
    KK = [];
    
    reg = cell(object.nDyn,1);
    
    for i=1:object.nDyn
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
    
    res = zeros(numel(Kred),object.nDyn);
    
    for i=1:object.nDyn
        centr = reg{i}.chebyCenter;
        if centr.r > 0
            res(:,i) = (Hred*centr.x <= Kred);
        else
            error('Empty dynamic found! Chech H and K matrix');
        end
    end
        fprintf(f,'int findDynamic(float inputVect[nIn_OBS], float state[nOut_OBS])\n{\n');
fprintf(f,'const float Hedge[%d][nx_OBS+np_OBS+nd_OBS] = {{',numel(Kred));
for i=1:numel(Kred)-1
    for j=1:size(Hred,2)-1
        fprintf(f,'%.10e,',Hred(i,j));
    end
    fprintf(f,'%.10e},\n{',Hred(i,end));
end
for j=1:size(Hred,2)-1
    fprintf(f,'%.10e,',Hred(end,j));
end
fprintf(f,'%.10e}};\n\n',Hred(end,end));

fprintf(f,'const float Kedge[%d] = {',numel(Kred));
for i=1:numel(Kred)-1
    fprintf(f,'%.10e,\n',Kred(i));
end
fprintf(f,'%.10e};\n\n',Kred(end));
fprintf(f,'int i,j,res;\n');
fprintf(f,'float tmp;\n');
fprintf(f,'res = 0;\n');
fprintf(f,'for(i=0;i<%d;i++)\n{\n',numel(Kred));
fprintf(f,'tmp = 0;\n');
fprintf(f,'for(j=0;j<nx_OBS;j++)\n');
fprintf(f,'    tmp += Hedge[i][j]*state[j];\n');
fprintf(f,'for(j=0;j<np_OBS;j++)\n');
fprintf(f,'    tmp += Hedge[i][j+nx_OBS]*inputVect[j+nu_OBS];\n');
fprintf(f,'for(j=0;j<nd_OBS;j++)\n');
fprintf(f,'    tmp += Hedge[i][j+nx_OBS+np_OBS]*state[j+nx_OBS];\n');
fprintf(f,'tmp -= Kedge[i];\n');
fprintf(f,'if (tmp <= 0)\n');
fprintf(f,'    res += (1<<i);\n');
fprintf(f,'}\n');    
intRes = zeros(1,object.nDyn);

for i=1:object.nDyn
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .H FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'observer.h'],'w');

fprintf(f,'#ifndef _OBSERVER_H_\n');
fprintf(f,'#define _OBSERVER_H_\n');
fprintf(f,'#define nIn_OBS %d\n',nIn);
fprintf(f,'#define nOut_OBS %d\n',nOut);
fprintf(f,'#define nDyn_OBS %d\n',object.nDyn);
fprintf(f,'#define nx_OBS %d\n',object.nx);
fprintf(f,'#define np_OBS %d\n',object.np);
fprintf(f,'#define nd_OBS %d\n',object.nd);
fprintf(f,'#define nu_OBS %d\n',object.nu);
fprintf(f,'void predict(float observerInput[nIn_OBS], float current_state[nOut_OBS], float predict_state[nOut_OBS],int dynamic);\n');
fprintf(f,'#endif\n');
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'observer.c'],'w');

fprintf(f,'#include "observer.h"\n');
[Aobs, Bobs, Cobs, Dobs, GxobsAll, Gyobs] = object.getMatrices();



fprintf(f,'const float Fpred[nDyn_OBS][nOut_OBS][nOut_OBS+nIn_OBS] = {{{');
for k=1:object.nDyn
    F = [Aobs{k} Bobs{k}];
    
    for i=1:nOut-1
        for j=1:nIn+nOut-1
            fprintf(f,'%.10e,',F(i,j));
        end
        fprintf(f,'%.10e},\n{',F(i,end));
    end
    for j=1:nIn+nOut-1
        fprintf(f,'%.10e,',F(end,j));
    end
    fprintf(f,'%.10e}}',F(end,end));
    
    if k ~= object.nDyn
        fprintf(f,',\n{{');
    else
        fprintf(f,'};\n\n');
    end
end


fprintf(f,'const float Gpred[nDyn_OBS][nOut_OBS] = {{');
for k=1:object.nDyn
    Gxobs = GxobsAll{k};
    for i=1:nOut-1
        fprintf(f,'%.10e,\n',Gxobs(i));
    end
    fprintf(f,'%.10e}',Gxobs(end));
    if k ~= object.nDyn
        fprintf(f,',\n{');
    else
        fprintf(f,'};\n\n');
    end
end



fprintf(f,'\n');
fprintf(f,'void predict(float observerInput[nIn_OBS], float current_state[nOut_OBS], float predict_state[nOut_OBS],int dynamic)\n{\n');
fprintf(f,'int i,j;\n');
fprintf(f,'for(i=0;i<nOut_OBS;i++)\n{\n');
fprintf(f,'\tpredict_state[i] = 0;\n');
fprintf(f,'\tfor(j=0;j<nOut_OBS;j++)\n');
fprintf(f,'\t\tpredict_state[i] += Fpred[dynamic][i][j]*current_state[j];\n');
fprintf(f,'\tfor(j=0;j<nIn_OBS;j++)\n');
fprintf(f,'\t\tpredict_state[i] += Fpred[dynamic][i][j+nOut_OBS]*observerInput[j];\n');
fprintf(f,'\tpredict_state[i] += Gpred[dynamic][i];\n');
fprintf(f,'}\n}\n');

fclose(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memory_size = numel([Aobs{1} Bobs{1} GxobsAll{1}])*object.nDyn;

if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    optO.memory_size = memory_size;
    optO.range = range;
    varargout{1} = optO;
end


% if needed generate interface
if circuit_parameters.generate_main
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
    for i = 1:nu
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:np
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nu),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nu),nbit,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:ny
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nu+np),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nu+np),nbit,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> y%d: [%s %s]\n',object.ynames{i},range.ymin(i),range.ymax(i),i,cirmin.bin,cirmax.bin);
    end
    
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    for i = 1:nx
        ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
        ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,ucirmin.bin,ucirmax.bin);
    end
    for i = 1:nd
        ucirmin = decimal2signed(circuit_parameters.outputRange.min(i+nx),nbitout,0);
        ucirmax = decimal2signed(circuit_parameters.outputRange.max(i+nx),nbitout,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,ucirmin.bin,ucirmax.bin);
    end
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])
    
    
    
end