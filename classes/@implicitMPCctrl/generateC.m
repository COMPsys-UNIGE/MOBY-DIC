function varargout = generateC(object,varargin)
% generateC  Generates C files for the circuit implementation of the MPC
%            controller on microcontroller
%
% generateC(OBJ)
% Generates the C files for the circuit implementation on microcontroller
% of implicit MPC controller with the default options.
% A report is also generated showing the main circuit features.
% The data are represented as float, therefore a 32bit floating point data
% representation is used.
%
% generateC(OBJ,OPTS)
% Generates the C files for the circuit implementation on microcontroller
% of the implicit MPC controller with custom options.
% OPTS is a structure with the following fields:
%  -  inputResolution: number of bits used to represent the function
%                      inputs. In general, it corresponds to the resolution
%                      of the analog-to-digital converters (ADC). Default
%                      value: 12.
%  -  inputRange: structure with fields min and max, indicating the minimum
%                 and maximum value of the circuit inputs. min and max must
%                 be integer values representable with the number of bits
%                 chosen in inputResolution. For example if inputResolution
%                 = 8 and inputRepresentation = 'unsigned', min and max
%                 must be integer numbers comprised betwen 0 and 255. In
%                 this case, if 0 and 255 are chosen as inputRange.min and
%                 inputRange.max, the maximum representation range is
%                 exploited. Default: maximum representable range.
%  - generateMain: flag indicating if the main file must be generated.
%                  If the controller is part of an embedded system, the main
%                  file must not be generated. Default value: true.
%  - folder: folder where C files are created. Default value: progressive folder.
%
% PERF = generateC(OBJ,...)
% If an output argument is provided, the report is not automatically opened
% and a structure PERF containing some information about the circuit
% implementation is returned. PERF is a structure with the following fields:
% TO DO

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
    circuit_parameters = [];
elseif nargin == 2
    circuit_parameters = varargin{1};
else
    error('Wrong input arguments');
end

circuit_parameters = MPCctrlCset(object,circuit_parameters);

% Number of ADMM iterations
maxIter = circuit_parameters.ADMMparameters.maxIter;

% Regularization parameter for ADMM
regPar = circuit_parameters.ADMMparameters.regPar;

% Number of variables
nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nxref = numel(object.getTrackingVariable);
trackVar = object.getTrackingVariable;
nu = object.getNumberOfInputs;
N = object.options.N;
Nu = object.options.Nu;

folder = circuit_parameters.folder;

%% Compute matrices for algorithm initializations

% LTI system matrices
augSys = object.getInformation.sys;
[A, ~, ~, ~, Ex, ~, Fx, ~, Gx, ~] = augSys.getMatrices;

if nxref ~= 0
    augSys = augSys.deltau();
    %     x = [x; uOld];
    % LTI extended system matrices
    [A, ~, ~, ~, Ex, ~, Fx, ~, Gx, ~] = augSys.getMatrices;
end

%% Compute matrices for ADMM algorithm
if ~object.options.tracking
    [QP1, QP2, ineqL, ineqR, eqL, eqR] = object.computeQP(zeros(nx,1), zeros(np,1), zeros(nd,1), object.options.ref, "ADMM");
else
    [QP1, QP2, ineqL, ineqR, eqL, eqR] = object.computeQP(zeros(nx,1), zeros(np,1), zeros(nd,1), zeros(nxref,1), "ADMM");
end
[M1, M2, M3, M4, v1, v2, v3] = object.computeADMMmatrices(QP1, QP2, ineqL, ineqR, eqL, eqR, regPar);


%% Constraints
% Get indices of states with soft constraints
[Hs,~] = object.constr.getAllConstraints('soft');
softIdx = [];
for i=1:size(Hs{1},2)
    if any(Hs{1}(:,i)~=0)
        softIdx = [softIdx, i];
    end
end
nSoft = length(softIdx);

% ----- TO DO: ADD SOFT CONSTRAINTS FOR ADMM -----
% [~,Kx] = object.constr.getStateConstraints;
% if nSoft == 0
%     range.xmin = -Kx{1}(1:nx);
%     range.xmax = Kx{1}(nx+1:2*nx);
% else
%     range.xmin = [];
%     range.xmax = [];
%     for i=1:2:2*(nx)
%         range.xmin = [range.xmin; -Kx{1}(i)];
%         range.xmax = [range.xmax; Kx{1}(i+1)];
%     end
% end
% 
% [~,Kp] = object.constr.getParameterConstraints;
% for i=1:np
%     range.pmin(i) = -Kp{1}(i);
%     range.pmax(i) = Kp{1}(np+i);
% end
% 
% [~,Kd] = object.constr.getUnmeasurableInputConstraints;
% for i=1:nd
%     range.dmin(i) = -Kd{1}(i);
%     range.dmax(i) = Kd{1}(nd+i);
% end
% 
% [~,Kxref] = object.constr.getReferenceConstraints;
% for i=1:nxref
%     range.xrefmin(i) = -Kxref{1}(i);
%     range.xrefmax(i) = Kxref{1}(nxref+i);
% end
% 
% [~,Ku] = object.constr.getInputConstraints;
% for i=1:nu
%     range.umin(i) = -Ku{1}(i);
%     range.umax(i) = Ku{1}(nu+i);
% end

range = circuit_parameters.range;

disp(['Destination folder for C files: ', folder]);
disp('Generating C files...');

if ~exist(folder,'dir')
    mkdir(folder);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    fprintf(f,'void scaleInput(int inADC[nX+nP+nD+nRef], float inscale[nX+nP+nD+nRef]);\n');
    fprintf(f,'void scaleOutput(int uADC[nU], float uscale[nU]);\n');
    fprintf(f,'\n\n');
    fprintf(f,'char controlGo = 0;\n');
    fprintf(f,'\n\n');
    fprintf(f,'void interruptTimerRoutine()\n{\n');
    fprintf(f,'controlGo = 1;\n');
    fprintf(f,'}\n\n');
    fprintf(f,'\n\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'\tint inADC[nX+nP+nD+nRef];//input vector\n');
    fprintf(f,'\tint inscale[nX+nP+nD+nRef];//scaled input vector\n');
    fprintf(f,'\tfloat x[nX];//state\n');

    if np ~= 0
        fprintf(f,'\tfloat p[nP];//parameter\n');
    end
    if nd ~= 0
        fprintf(f,'\tfloat d[nD];//unmeasurable input\n');
    end
    if nxref ~= 0
        fprintf(f,'\tfloat ref[nRef];//reference state\n');
    end

    fprintf(f,'\tfloat u[nU];//vector containing control value \n\n');
    fprintf(f,'\tint uADC[nU]\n');

    fprintf(f,'\tconfigADC();//config ADC \n');
    fprintf(f,'\tconfigDAC();//config DAC  \n\n');
    fprintf(f,'\tconfigTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',object.getSamplingTime);
    fprintf(f,'\twhile(1)\n\t{\n');
    fprintf(f,'\t\twhile(!controlGo);//wait for interrupt activation\n');
    fprintf(f,'\t\tcontrolGo = 0;\n');
    fprintf(f,'\t\treadInputADC(inADC);\n');

    fprintf(f,'\t\tscaleInput(inADC,inscale);\n');
    for i=1:nx
        fprintf(f,'\t\tx[%d] = inscale[%d];\n',i-1,i-1);
    end
    for i=nx+1:nx+np
        fprintf(f,'\t\tp[%d] = inscale[%d];\n',i-nx-1,i-1);
    end
    for i=nx+np+1:nx+np+nd
        fprintf(f,'\t\td[%d] = inscale[%d];\n',i-nx-np-1,i-1);
    end
    for i=nx+np+nd+1:nx+np+nd+nxref
        fprintf(f,'\t\tref[%d] = inscale[%d];\n',i-nx-np-nd-1,i-1);
    end

    if nxref ~= 0
        if np ~= 0 && nd ~= 0
            fprintf(f,'\t\tcontrol(x,u,p,d,ref);\n');
        elseif np ~= 0
            fprintf(f,'\t\tcontrol(x,u,p,ref);\n');
        elseif nd ~= 0
            fprintf(f,'\t\tcontrol(x,u,d,ref);\n');
        else
            fprintf(f,'\t\tcontrol(x,u,ref);\n');
        end
    else
        if np ~= 0 && nd ~= 0
            fprintf(f,'\t\tcontrol(x,u,p,d);\n');
        elseif np ~= 0
            fprintf(f,'\t\tcontrol(x,u,p);\n');
        elseif nd ~= 0
            fprintf(f,'\t\tcontrol(x,u,d);\n');
        else
            fprintf(f,'\t\tcontrol(x,u);\n');
        end
    end

    fprintf(f,'\t\tscaleOutput(uADC,u);\n');
    fprintf(f,'\t\tsetOutputDAC(uADC);\n\t}\n');
    fprintf(f,'}\n\n');

    fprintf(f,'void scaleInput(int inADC[nX+nP+nD+nRef], float inscale[nX+nP+nD+nRef])\n{\n');
    for i=1:nx
        deltaReal = range.xmax(i)-range.xmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.xmax(i)+range.xmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        fprintf(f,'\tinscale[%d] = (%.10e)*inADC[%d]+(%.10e);\n',i-1,deltaReal/deltaADC,i-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:np
        deltaReal = range.pmax(i)-range.pmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.pmax(i)+range.pmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        fprintf(f,'\tinscale[%d] = (%.10e)*inADC[%d]+(%.10e);\n',i+nx-1,deltaReal/deltaADC,i+nx-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:nd
        deltaReal = range.dmax(i)-range.dmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.dmax(i)+range.dmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        fprintf(f,'\tinscale[%d] = (%.10e)*inADC[%d]+(%.10e);\n',i+nx+np-1,deltaReal/deltaADC,i+nx+np-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    for i=1:nxref
        deltaReal = range.xrefmax(i)-range.xrefmin(i);
        deltaADC = circuit_parameters.inputRange.max(i)-circuit_parameters.inputRange.min(i);
        meanReal = (range.xrefmax(i)+range.xrefmin(i))/2;
        meanADC = (circuit_parameters.inputRange.max(i)+circuit_parameters.inputRange.min(i))/2;
        fprintf(f,'\tinscale[%d] = (%.10e)*inADC[%d]+(%.10e);\n',i+nx+np+nd-1,deltaReal/deltaADC,i+nx+np+nd-1,(meanReal-deltaReal/deltaADC*meanADC));
    end
    fprintf(f,'}\n');

    fprintf(f,'void scaleOutput(int uADC[nU], float uscale[nU])\n{\n');
    for i=1:nu
        deltaReal =  range.umax(i)-range.umin(i);
        deltaADC = circuit_parameters.outputRange.max(i)-circuit_parameters.outputRange.min(i);
        meanReal = (range.umax(i)+range.umin(i))/2;
        meanADC = (circuit_parameters.outputRange.max(i)+circuit_parameters.outputRange.min(i))/2;
        fprintf(f,'\tuADC[%d] = (%.10e)*uscale[%d]+(%.10e);\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    fclose(f);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% GENERATE CONTROLLER.H FILE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----- TO DO: ADD SOFT CONSTRAINTS FOR ADMM -----
% Get number of states with soft constraints
% [Hs,~] = object.constr.getAllConstraints('soft');
% softIdx = [];
% for i=1:size(Hs{1},2)
%     if any(Hs{1}(:,i)~=0)
%         softIdx = [softIdx, i];
%     end
% end
% nSoft = length(softIdx);
% nSoft = 0;

% Open the header file
f = fopen([folder,'controller.h'],'w');

fprintf(f,'#ifndef _ADMM_H_\n');
fprintf(f,'#define _ADMM_H_\n\n');

fprintf(f,'// Number of system states (x_k)\n');
fprintf(f,'#define nX %d\n\n',nx);

fprintf(f,'// Number of system inputs (u_k)\n');
fprintf(f,'#define nU %d\n\n',nu);

fprintf(f,'// Number of system parameters (p_k)\n');
fprintf(f,'#define nP %d\n\n',np);

fprintf(f,'// Number of system unmeasurable inputs (d_k)\n');
fprintf(f,'#define nD %d\n\n',nd);

fprintf(f,'// Prediction horizon\n');
fprintf(f,'#define N %d\n\n',N);

fprintf(f,'// Control horizon\n');
fprintf(f,'#define Nu %d\n\n',Nu);

if nxref == 0
    nDim = Nu*nu + (N+1)*nx;
    fprintf(f,'// Number of controller inputs [x_k]\n');
    fprintf(f,'#define nX_CTRL %d\n\n',nx);
else
    fprintf(f,'// Number of reference inputs\n');
    fprintf(f,'#define nRef %d\n\n',nxref);
    nDim = Nu*nu + (N+1)*(nx+nu);
    fprintf(f,'// Number of controller inputs [x_k; u_(k-1)]\n');
    fprintf(f,'#define nX_CTRL %d\n\n',nx+nu);
end

fprintf(f,'// Number of optimization variables\n');
fprintf(f,'#define nDim_CTRL %d\n\n', nDim);

nCon = numel(v3);
fprintf(f,'// Number of constraints\n');
fprintf(f,'#define nCon %d\n\n', nCon);

fprintf(f,'// Regularization parameter\n');
fprintf(f,'#define reg_par %d\n\n', regPar);

fprintf(f,'// Maximum number of iterations\n');
fprintf(f,'#define max_iter %d\n\n', maxIter);

if nxref ~= 0
    fprintf(f,'// Index of the input to be tracked\n');
    fprintf(f,'static const float ref_idx[nRef] = {');
    for i=1:nxref-1
        fprintf(f,'%d, ', fix(trackVar(i)-1));
    end
    fprintf(f,'%d};\n\n', fix(trackVar(nxref)-1));
end

fprintf(f,'// Precomputed matrix for ADMM initialization\n');
if nxref ~= 0
    Apwr = [];
    for i=1:N
       Apwr = [Apwr; A^i];
    end
    fprintf(f,'static const float A_pwr[N*nX_CTRL][nX_CTRL] = {\n');
    for i=1:N*(nx+nu)-1
        fprintf(f,'{%f',Apwr(i,1));
        for j=2:nx+nu
            fprintf(f,', ');
            fprintf(f,'%f',Apwr(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',Apwr(N*(nx+nu),1));
    for j=2:nx+nu
        fprintf(f,', ');
        fprintf(f,'%f',Apwr(N*(nx+nu),j));
    end
else
    Apwr = [];
    for i=1:N
       Apwr = [Apwr; A^i];
    end
    fprintf(f,'static const float A_pwr[N*nX_CTRL][nX_CTRL] = {\n');
    for i=1:N*nx-1
        fprintf(f,'{%f',Apwr(i,1));
        for j=2:nx
            fprintf(f,', ');
            fprintf(f,'%f',Apwr(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',Apwr(N*(nx),1));
    for j=2:nx
        fprintf(f,', ');
        fprintf(f,'%f',Apwr(N*(nx),j));
    end
end
fprintf(f,'} };\n\n');

if np ~= 0
    fprintf(f,'// Precomputed matrix for ADMM initialization\n');
    fprintf(f,'static const float E[nX_CTRL][nP] = {\n');
    if nxref ~= 0
        for i=1:nx+nu-1
          fprintf(f,'{%f', Ex(i,1));
          for j=2:np
              fprintf(f,', ');
              fprintf(f,'%f', Ex(i,j));
          end
          fprintf(f,'},\n');
        end
        fprintf(f,'{%f', Ex(nx+nu,1));
        for j=2:np
            fprintf(f,', ');
            fprintf(f,'%f', Ex(nx+nu,j));
        end
    else
        for i=1:nx-1
          fprintf(f,'{%f', Ex(i,1));
          for j=2:np
              fprintf(f,', ');
              fprintf(f,'%f', Ex(i,j));
          end
          fprintf(f,'},\n');
        end
        fprintf(f,'{%f', Ex(nx,1));
        for j=2:np
            fprintf(f,', ');
            fprintf(f,'%f', Ex(nx,j));
        end
    end
    fprintf(f,'} };\n\n');
end

if nd ~= 0
    fprintf(f,'// Precomputed matrix for ADMM initialization\n');
    fprintf(f,'static const float F[nX_CTRL][nD] = {\n');
    if nxref ~= 0
        for i=1:nx+nu-1
          fprintf(f,'{%f', Fx(i,1));
          for j=2:nd
              fprintf(f,', ');
              fprintf(f,'%f', Fx(i,j));
          end
          fprintf(f,'},\n');
        end
        fprintf(f,'{%f', Fx(nx+nu,1));
        for j=2:np
            fprintf(f,', ');
            fprintf(f,'%f', Fx(nx+nu,j));
        end
    else
        for i=1:nx-1
          fprintf(f,'{%f', Fx(i,1));
          for j=2:nd
              fprintf(f,', ');
              fprintf(f,'%f', Fx(i,j));
          end
          fprintf(f,'},\n');
        end
        fprintf(f,'{%f', Fx(nx,1));
        for j=2:np
            fprintf(f,', ');
            fprintf(f,'%f', Fx(nx,j));
        end
    end
    fprintf(f,'} };\n\n');
end

fprintf(f,'// Precomputed matrix for ADMM initialization\n');
if nxref ~= 0
    fprintf(f,'static const float G[nX_CTRL] = {');
    for i=1:nx+nu-1
        fprintf(f,'%f, ', Gx(i));
    end
    fprintf(f,'%f};\n\n', Gx(nx+nu));
else
    fprintf(f,'static const float G[nX_CTRL] = {');
    for i=1:nx-1
        fprintf(f,'%f, ', Gx(i));
    end
    fprintf(f,'%f};\n\n', Gx(nx));
end

fprintf(f,'// Precomputed matrix for ADMM (it depends on A and B system matrices)\n');
fprintf(f,'static const float M1[nDim_CTRL][nDim_CTRL] = {\n');
for i=1:nDim-1
    fprintf(f,'{%f',M1(i,1));
    for j=2:nDim
        fprintf(f,', ');
        fprintf(f,'%f',M1(i,j));
    end
    fprintf(f,'},\n');
end
fprintf(f,'{%f',M1(nDim,1));
for j=2:nDim
    fprintf(f,', ');
    fprintf(f,'%f',M1(nDim,j));
end
fprintf(f,'} };\n\n');

fprintf(f,'// Precomputed matrix for ADMM\n');
fprintf(f,'static const float M2[nDim_CTRL][nCon + N*nX_CTRL] = {\n');
if nxref ~= 0
    for i=1:nDim-1
        fprintf(f,'{%f',M2(i,1));
        for j=2:nCon + N*(nx+nu)
            fprintf(f,', ');
            fprintf(f,'%f',M2(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',M2(nDim,1));
    for j=2:nCon + N*(nx+nu)
        fprintf(f,', ');
        fprintf(f,'%f',M2(nDim,j));
    end
else
    for i=1:nDim-1
        fprintf(f,'{%f',M2(i,1));
        for j=2:nCon + N*nx
            fprintf(f,', ');
            fprintf(f,'%f',M2(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',M2(nDim,1));
    for j=2:nCon + N*nx
        fprintf(f,', ');
        fprintf(f,'%f',M2(nDim,j));
    end
end
fprintf(f,'} };\n\n');

fprintf(f,'// Precomputed matrix for ADMM\n');
fprintf(f,'static const float M3[nCon][nDim_CTRL] = {\n');
for i=1:nCon-1
    fprintf(f,'{%f',M3(i,1));
    for j=2:nDim
        fprintf(f,', ');
        fprintf(f,'%f',M3(i,j));
    end
    fprintf(f,'},\n');
end
fprintf(f,'{%f',M3(nCon,1));
for j=2:nDim
    fprintf(f,', ');
    fprintf(f,'%f',M3(nCon,j));
end
fprintf(f,'} };\n\n');

fprintf(f,'// Precomputed matrix for ADMM\n');
fprintf(f,'static const float M4[nCon + N*nX_CTRL][nDim_CTRL] = {\n');
if nxref ~= 0
    for i=1:nCon + N*(nx+nu) -1
        fprintf(f,'{%f',M4(i,1));
        for j=2:nDim
            fprintf(f,', ');
            fprintf(f,'%f',M4(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',M4(nCon + N*(nx+nu),1));
    for j=2:nDim
        fprintf(f,', ');
        fprintf(f,'%f',M4(nCon + N*(nx+nu),j));
    end
else
    for i=1:nCon + N*nx -1
        fprintf(f,'{%f',M4(i,1));
        for j=2:nDim
            fprintf(f,', ');
            fprintf(f,'%f',M4(i,j));
        end
        fprintf(f,'},\n');
    end
    fprintf(f,'{%f',M4(nCon + N*nx,1));
    for j=2:nDim
        fprintf(f,', ');
        fprintf(f,'%f',M4(nCon + N*nx,j));
    end
end
fprintf(f,'} };\n\n');

fprintf(f,'// Precomputed vector for ADMM\n');
fprintf(f,'static const float v3[nCon] = {');
for i=1:nCon - 1
    fprintf(f,'%f, ', v3(i));
end
fprintf(f,'%f};\n\n', v3(nCon));

if nxref ~= 0
    fprintf(f,'// Weight matrix for the state vector in the optimization\n');
    fprintf(f,'// problem (useful when the reference value changes)\n');
    Q = zeros(nx+nu);
    Q(1:nx,1:nx) = object.options.Q;
    fprintf(f,'static const float Q[nX_CTRL][nX_CTRL] = {');
    for i=1:nx+nu-1
        fprintf(f,'{%f',Q(i,1));
        for j=2:nx+nu
            fprintf(f,', ');
            fprintf(f,'%f',Q(i,j));
        end
    fprintf(f,'},\n');
    end
    fprintf(f,'{%f',Q(nx+nu,1));
    for j=2:nx+nu
        fprintf(f,', ');
        fprintf(f,'%f',Q(nx+nu,j));
    end
    fprintf(f,'} };\n\n');

    fprintf(f,'// Weight matrix for the state vector in the optimization\n');
    fprintf(f,'// problem (useful when the reference value changes)\n');
    P = zeros(nx+nu);
    P(1:nx,1:nx) = object.options.P;
    fprintf(f,'static const float P[nX_CTRL][nX_CTRL] = {');
    for i=1:nx+nu-1
        fprintf(f,'{%f',P(i,1));
        for j=2:nx+nu
            fprintf(f,', ');
            fprintf(f,'%f',P(i,j));
        end
    fprintf(f,'},\n');
    end
    fprintf(f,'{%f',P(nx+nu,1));
    for j=2:nx+nu
        fprintf(f,', ');
        fprintf(f,'%f',P(nx+nu,j));
    end
    fprintf(f,'} };\n\n');
else
    fprintf(f,'static const float q[nDim_CTRL] = {');
    for i=1:nDim-1
        fprintf(f,'%f, ',QP2(i));
    end
    fprintf(f,'%f};\n\n',QP2(nDim));
end

if nxref ~= 0
    fprintf(f,'// Default control\n');
    fprintf(f,'static const float default_u[nU] = {');
    for i=1:nu - 1
%         fprintf(f,'%f, ',do(i));
        fprintf(f,'%f, ',0);
    end
%     fprintf(f,'%f};\n\n',do(nu));
    fprintf(f,'%f};\n\n',0);
end

if nxref ~= 0
    fprintf(f,'void augmentState(float x_reg[nX], float ref_reg[nX], float x_aug[nX_CTRL], float ref_aug[nX_CTRL], float u_old[nU]);\n');
    fprintf(f,'void extractU(float z[nDim_CTRL], float u_reg[nU], float u_old[nU]);\n');
    if np ~= 0 && nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP], float d_in[nD], float ref_in[nRef]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float p[nP], float d[nD], float ref[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    elseif np ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP], float ref_in[nRef]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float p[nP], float ref[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    elseif nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float d_in[nD], float ref_in[nRef]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float d[nD], float ref[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    else
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float ref_in[nRef]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float ref[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    end
else
    if np ~= 0 && nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP], float d_in[nD]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float p[nP], float d[nD], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    elseif np ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float p[nP], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    elseif nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float d_in[nD]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float d[nD], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    else
        fprintf(f,'void control(float x_in[nX], float u_opt[nU]);\n\n');
        fprintf(f,'void admmInit(float x[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL]);\n');
    end
end
fprintf(f,'void admm(float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL], float z[nDim_CTRL]);\n');

fprintf(f,'#endif\n');

% Close the header file
fclose(f);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% GENERATE CONTROLLER.C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open the main file
f = fopen([folder,'controller.c'],'w');

fprintf(f,'#include "controller.h"\n\n');

if nxref ~= 0
    if np ~= 0 && nd ~= 0
        fprintf(f,'void control(float x_in[nX_CTRL], float u_opt[nU], float p_in[nP], float d_in[nD], float ref_in[nRef])\n{\n');
    elseif np ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP], float ref_in[nRef])\n{\n');
    elseif nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float d_in[nD], float ref_in[nRef])\n{\n');
    else
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float ref_in[nRef])\n{\n');
    end
else
    if np ~= 0 && nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP], float d_in[nD])\n{\n');
    elseif np ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float p_in[nP])\n{\n');
    elseif nd ~= 0
        fprintf(f,'void control(float x_in[nX], float u_opt[nU], float d_in[nD])\n{\n');
    else
        fprintf(f,'void control(float x_in[nX], float u_opt[nU])\n{\n');
    end
end

fprintf(f,'float x_reg[nX_CTRL];\n\n');
fprintf(f,'float u_reg[nU];\n\n');
if nxref ~= 0
    fprintf(f,'float x_aug[nX_CTRL];\n\n');
    fprintf(f,'float ref_reg[nX];\n\n');
    fprintf(f,'float ref_aug[nX_CTRL];\n\n');
    fprintf(f,'static float u_old[nU];\n\n');
    fprintf(f,'float delta_u[nU];\n\n');
    fprintf(f,'static bool firstCall = true;\n\n');
    fprintf(f,'static float q[nDim_CTRL];\n\n');
end
if np ~= 0
    fprintf(f,'float p_reg[nP];\n\n');
end
if nd ~= 0
    fprintf(f,'float d_reg[nD];\n\n');
end

fprintf(f,'static float v1[nDim_CTRL];\n\n');
fprintf(f,'static float v2[nCon + N*nX_CTRL];\n\n');

fprintf(f,'float z[nDim_CTRL];\n\n');

fprintf(f,'for (int i = 0; i < nX; i++)\n{\n');
fprintf(f,'\tx_reg[i] = x_in[i];\n');
if nxref ~= 0
    fprintf(f,'\tref_reg[i] = 0;\n}\n\n');
    fprintf(f,'for (int i = 0; i < nRef; i++)\n{\n');
    fprintf(f,'\tref_reg[ref_idx[i]] = ref_in[i];\n}\n\n');
else
    fprintf(f,'}\n\n');
end
if np ~= 0
    fprintf(f,'for (int i = 0; i < nP; i++)\n{\n');
    fprintf(f,'\tp_reg[i] = p_in[i];\n}\n\n');
end
if nd ~= 0
    fprintf(f,'for (int i = 0; i < nD; i++)\n{\n');
    fprintf(f,'\td_reg[i] = d_in[i];\n}\n\n');
end

if nxref ~= 0
    fprintf(f,'if (firstCall)\n{\n');
    fprintf(f,'\tfor (int i = 0; i < nU; i++)\n{\n');
    fprintf(f,'\tu_old[i] = default_u[i];\n}\n');
    fprintf(f,'\tfirstCall = false;\n}\n\n');
    fprintf(f,'augmentState(x_reg, ref_reg, x_aug, ref_aug, u_old);\n\n');

    if np ~= 0 && nd ~= 0
        fprintf(f,'admmInit(x_aug, p_reg, d_reg, ref_aug, v1, v2);\n\n');
    elseif np ~= 0
        fprintf(f,'admmInit(x_aug, p_reg, ref_aug, v1, v2);\n\n');
    elseif nd ~= 0
        fprintf(f,'admmInit(x_aug, d_reg, ref_aug, v1, v2);\n\n');
    else
        fprintf(f,'admmInit(x_aug, ref_aug, v1, v2);\n\n');
    end
    fprintf(f,'admm(v1, v2, z);\n\n');
    
    fprintf(f,'extractU(z, u_reg, u_old);\n\n');
else
    if np ~= 0 && nd ~= 0
        fprintf(f,'admmInit(x_reg, p_reg, d_reg, v1, v2);\n\n');
    elseif np ~= 0
        fprintf(f,'admmInit(x_reg, p_reg, v1, v2);\n\n');
    elseif nd ~= 0
        fprintf(f,'admmInit(x_reg, d_reg, v1, v2);\n\n');
    else
        fprintf(f,'admmInit(x_reg, v1, v2);\n\n');
    end
    fprintf(f,'admm(v1, v2, z);\n\n');
    
    fprintf(f,'for (int i = 0; i < nU; i++)\n{\n');
    fprintf(f,'\t#pragma HLS UNROLL\n');
    fprintf(f,'\tu_reg[i] = z[i];\n}\n\n');
end

fprintf(f,'for (int i = 0; i < nU; i++)\n{\n');
fprintf(f,'\t#pragma HLS UNROLL\n');
fprintf(f,'\tu_opt[i] = u_reg[i];\n}\n\n');

fprintf(f,'}\n');

% Close the main file
fclose(f);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% COPY C FILES FOR CONTROLLER AND ADMM %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy the appropriate C files
copyfile([getcpath,'/implicitMPCctrl/admmInit.c'],[folder,'admmInit.c']);
copyfile([getcpath,'/implicitMPCctrl/admm.c'],[folder,'admm.c']);

if nxref ~= 0
    if np ~= 0 && nd ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd p[nP], fxd d[nD], fxd ref[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    elseif np ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd p[nP], fxd ref[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    elseif nd ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd d[nD], fxd ref[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    else
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd ref[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    end
else
    if np ~= 0 && nd ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd p[nP], fxd d[nD], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    elseif np ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd p[nP], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    elseif nd ~= 0
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd d[nD], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    else
        S = fileread([folder,'admmInit.c']);
        S = strrep(S,'void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])','void admmInit(fxd x[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + N*nX_CTRL])');
        f = fopen([folder,'admmInit.c'],'w');
        fwrite(f, S);
        fclose(f);
    end
end

if np ~= 0
    S = fileread([folder,'admmInit.c']);
    S = strrep(S,'/*p','');
    S = strrep(S,'p*/','');
    f = fopen([folder,'admmInit.c'],'w');
    fwrite(f, S);
    fclose(f);
end

if nd ~= 0
    S = fileread([folder,'admmInit.c']);
    S = strrep(S,'/*d','');
    S = strrep(S,'d*/','');
    f = fopen([folder,'admmInit.c'],'w');
    fwrite(f, S);
    fclose(f);
end

if nxref ~= 0
    S = fileread([folder,'admmInit.c']);
    S = strrep(S,'/*ref','');
    S = strrep(S,'ref*/','');
    f = fopen([folder,'admmInit.c'],'w');
    fwrite(f, S);
    fclose(f);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional function output
if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    optO.range = range;
    varargout{1} = optO;
end

% If an interface was generated, generate a report
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
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:np
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx),nbit,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:nd
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np),nbit,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,cirmin.bin,cirmax.bin);
    end
    for i = 1:nxref
        track = object.getTrackingVariable;
        cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np+nd),nbit,0);
        cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np+nd),nbit,0);
        fprintf(fout,'\t\t%s_ref: [%.10e %.10e] --> xref%d: [%s %s]\n',object.xnames{track(i)},range.xrefmin(i),range.xrefmax(i),i,cirmin.bin,cirmax.bin);
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
    fprintf(fout,'\t - Range (model --> circuit):\n');
    for i = 1:nu
        ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
        ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
        fprintf(fout,'\t\t%s: [%.10e %.10e] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
    end

    fclose(fout);
    edit([folder ,'C_report.log'])

end
