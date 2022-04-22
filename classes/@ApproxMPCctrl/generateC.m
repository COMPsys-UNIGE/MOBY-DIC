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

circuit_parameters = object.ApproxMPCctrlCset(circuit_parameters);

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nxref = numel(object.getTrackingVariable);

nu = object.getNumberOfInputs;
nDim = nx+nd+np+nxref;
% Set default value for idx
idx = 1:nu;

% Get PWAG or PWAS function
ctrl = object;
pwaFun = ctrl.getFunction;

folder = circuit_parameters.folder;
subOpt.generate_main = false;
subOpt.folder = folder;
funReport = pwaFun.generateC(idx,subOpt);



range.xmin = funReport.range.xmin(1:nx);
range.xmax = funReport.range.xmax(1:nx);

range.pmin = funReport.range.xmin(nx+1:nx+np);
range.pmax = funReport.range.xmax(nx+1:nx+np);

range.dmin = funReport.range.xmin(nx+np+1:nx+np+nd);
range.dmax = funReport.range.xmax(nx+np+1:nx+np+nd);

range.xrefmin = funReport.range.xmin(nx+np+nd+1:end);
range.xrefmax = funReport.range.xmax(nx+np+nd+1:end);

range.umin = funReport.range.umin;
range.umax = funReport.range.umax;

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
    fprintf(f,'#include "controller.h"\n');
    fprintf(f,'void scaleInput(int xADC[nDim_CTRL], float xscale[nDim_CTRL]);\n');
    fprintf(f,'void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL]);\n');
    fprintf(f,'\n\n');
    fprintf(f,'char controlGo = 0;\n');
    fprintf(f,'\n\n');
    fprintf(f,'void interruptTimerRoutine()\n{\n');
    fprintf(f,'controlGo = 1;\n');
    fprintf(f,'}\n\n');
    fprintf(f,'\n\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'int xADC[nDim_CTRL]\n');
    fprintf(f,'float x[nDim_CTRL];//input vector\n');
    fprintf(f,'float u[nY_CTRL];//vector containing control value \n\n');
    fprintf(f,'int uADC[nY_CTRL]\n');
    fprintf(f,'configADC();//config ADC \n');
    fprintf(f,'configDAC();//config DAC  \n\n');
    fprintf(f,'configTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',object.getSamplingTime);
    fprintf(f,'while(1)\n{\n');
    fprintf(f,'while(!controlGo);//wait for interrupt activation\n');
    fprintf(f,'controlGo = 0;\n');
    fprintf(f,'readInputADC(xADC);\n');
    fprintf(f,'scaleInput(xADC,x);\n');
    fprintf(f,'/* x = [x1,x2,...,xn,p1,..,pq,d1,..,dz] */\n\n');
    fprintf(f,'control(x,u);\n');
    fprintf(f,'scaleOutput(uADC,u);\n');
    fprintf(f,'setOutputDAC(uADC);\n}\n');
    fprintf(f,'}\n\n');
    fprintf(f,'void scaleInput(int xADC[nDim_CTRL], float xscale[nDim_CTRL])\n{\n');
    fprintf(f,'int i;\n');
    fprintf(f,'for(i=0;i<nDim_CTRL;i++)\n');
    for i=1:nDim
        deltaReal =  funReport.range.xmax(i)-funReport.range.xmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (funReport.range.xmax(i)+funReport.range.xmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        
        fprintf(f,'xscale[i] = %.10e*xADC[i]+%.10e;\n',deltaReal/deltaADC,(meanReal-deltaReal/deltaADC*meanADC));
    end
    fprintf(f,'}\n');
    
    fprintf(f,'void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL])\n{\n');
    fprintf(f,'int i;\n');
    fprintf(f,'for(i=0;i<nY_CTRL;i++)\n');
    for i=1:nu
        deltaReal =  funReport.range.umax(i)-funReport.range.umin(i);
        deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
        meanReal = (funReport.range.umax(i)+funReport.range.umin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;
        
        fprintf(f,'uADC[i] = %.10e*uscale[i]+%.10e;\n',deltaADC/deltaReal,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .H FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'controller.h'],'w');

fprintf(f,'#ifndef _CONTROLLER_H_\n');
fprintf(f,'#define _CONTROLLER_H_\n');
fprintf(f,'#define nDim_CTRL %d\n',nDim);
fprintf(f,'#define nY_CTRL %d\n',nu);
fprintf(f,'void control(float *x, float *u);\n');
fprintf(f,'#endif\n');
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'controller.c'],'w');

fprintf(f,'#include "controller.h"\n');
fprintf(f,'#include "pwas_fun.h"\n\n');
fprintf(f,'void control(float *x, float *u)\n{\n');
fprintf(f,'calculatePWAS(x,u);\n');
fprintf(f,'}\n');
fclose(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memory_size = funReport.memory_size;

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
    for i = 1:nx
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
        fprintf(fout,'\t\t%s: [%f %f] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:np
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx),nbit,0);
        fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:nd
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np),nbit,0);
        fprintf(fout,'\t\t%s: [%f %f] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:nxref
        track = object.getTrackingVariable;
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np+nd),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np+nd),nbit,0);
        fprintf(fout,'\t\t%s_ref: [%f %f] --> xref%d: [%s %s]\n',object.xnames{track(i)},range.xrefmin(i),range.xrefmax(i),i,cirmin.bin,cirmax.bin);
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    for i = 1:nu
        ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
        ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
    end
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])
    
    
    
end