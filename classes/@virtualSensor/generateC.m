%% generateC
% Generates the C files describing the digital circuit implementing
% the MPCctrl
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

circuit_parameters = vsCset(object,circuit_parameters);

nu = object.getNumberOfInputs();
ny = object.getNumberOfMeasurableOutputs();
nz = object.getNumberOfUnmeasurableOutputs();

mu = object.getInputTimeWindow;
my = object.getOutputTimeWindow;
mz = object.getAutoregressiveTimeWindow;

domain = object.getDomain;
reducedComplexity = object.isReducedComplexity();

fun = object.getFunction();
npwas = numel(fun);

folder = circuit_parameters.folder;

% Current folder
thisfolder = pwd;

for i = 1:npwas
    curfolder = [folder,'pwas',num2str(i)];
    subOpt.folder = curfolder;
    subOpt.generate_main = false;
    
    funReport(i) = fun(i).generateC(1:nz,subOpt);
    
    % Move to folder containing the VHDL files of ii-th function
    cd(curfolder)
    
    % Rename file names and entity names by adding a progressive number to
    % the name itself; this allows to avoid name collisions in the circuit
    % containing all pwas functions
    files = dir;
    files(1:2) = [];
    nfiles = numel(files);
    filenamesc = cell(nfiles,1);
    filenamesh = cell(nfiles,1);
    
    kc = 1;
    kh = 1;
    for j = 1:nfiles
        tmpname = files(j).name;
        if strfind(tmpname,'.c')
            filenamesc{kc} = files(j).name(1:end-2);
            kc = kc+1;
        end
        if strfind(tmpname,'.h')
            filenamesh{kh} = files(j).name(1:end-2);
            kh = kh+1;
        end
    end
    filenamesc = filenamesc(1:kc-1);
    filenamesh = filenamesc(1:kh-1);
    nfilesc = numel(filenamesc);
    nfilesh = numel(filenamesh);
    
    namestochange{1} = 'pwas_fun';
    namestochange{2} = 'nDim';
    namestochange{3} = 'nY';
    namestochange{4} = 'nWeight';
    namestochange{5} = 'getMu';
    namestochange{6} = 'Sorting';
    namestochange{7} = 'Transform';
    namestochange{8} = 'divideZ';
    namestochange{9} = 'getSimplexVertices';
    namestochange{10} = 'calculatePWAS';
    namestochange{11} = 'uFunction';
    namestochange{12} = 'Npartition';
    namestochange{13} = 'mm';
    namestochange{14} = 'qq';
    namestochange{15} = 'weigthVector';
    namestochange{16} = '_PWAS_FUN_H';
      
    nnames = numel(namestochange);
    
    for j = 1:nfilesc
        fin = fopen([filenamesc{j},'.c'],'r');
        if fin ~= -1
            fout = fopen(['../',filenamesc{j},'_',num2str(i),'.c'],'w');
            
            while ~feof(fin)
                s = fgetl(fin);
                for k = 1:nnames
                    s = strrep(s, [namestochange{k}], [namestochange{k},'_',num2str(i)]);
                end
                fprintf(fout,'%s\n',s);
            end
            fclose(fin);
            fclose(fout);
        end
    end
    
    for j = 1:nfilesh
        fin = fopen([filenamesh{j},'.h'],'r');
        if fin ~= -1
            fout = fopen(['../',filenamesh{j},'_',num2str(i),'.h'],'w');
            
            while ~feof(fin)
                s = fgetl(fin);
                for k = 1:nnames
                    s = strrep(s, [namestochange{k}], [namestochange{k},'_',num2str(i)]);
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
    fprintf(f,'   readInputADC(float *y) - in this function the program must\n');
    fprintf(f,'         read value of all system output\n');
    fprintf(f,'   setOutputDAC() - in this function the program must write\n');
    fprintf(f,'         to DAC controls value\n');
    fprintf(f,'   It is also necessary to call function interruptTimerRoutine()\n');
    fprintf(f,'   into the timer interrupt routine.                            */\n');
    fprintf(f,'#include <stdio.h>\n');
    fprintf(f,'#include "virtualSensor.h"\n\n');
    fprintf(f,'void scaleInput(int uADC[nu_VS], float uscale[nu_VS]);\n');
    fprintf(f,'void scaleInput(int yADC[ny_VS], float yscale[ny_VS]);\n');
    fprintf(f,'void scaleOutput(int zADC[nz_VS], float zscale[nz_VS]);\n');
    fprintf(f,'\n\n');
    fprintf(f,'char controlGo = 0;\n');
    fprintf(f,'\n\n');
    fprintf(f,'void interruptTimerRoutine()\n{\n');
    fprintf(f,'\tcontrolGo = 1;\n');
    fprintf(f,'}\n\n');
    fprintf(f,'\n\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'\tint i,j;\n');
    fprintf(f,'\tint uyADC[nu_VS+ny_VS];\n');
    fprintf(f,'\tint zADC[nz_VS];\n');
    fprintf(f,'\tfloat uy[nu_VS+ny_VS];\n');
    fprintf(f,'\tfloat u[nu_VS]; //input vector\n');
    fprintf(f,'\tfloat y[ny_VS]; //measurable output vector\n');
    fprintf(f,'\tfloat z[nz_VS]; //unmeasurable output vector\n\n');
    fprintf(f,'\tfloat buf_u[nu_VS][mu_VS]; //buffer for inputs\n');
    fprintf(f,'\tfloat buf_y[ny_VS][my_VS]; //buffer for measurable outputs\n');
    fprintf(f,'\tfloat buf_z[nz_VS][mz_VS]; //buffer for unmeasurable outputs\n\n');
    
    fprintf(f,'\tint count_buf_u = 0;\n');
    fprintf(f,'\tint count_buf_y = 0;\n');
    fprintf(f,'\tint count_buf_z = 0;\n');
    fprintf(f,'\tint full_buf = 0;\n\n');
    
    fprintf(f,'\tconfigADC();//config ADC \n');
    fprintf(f,'\tconfigDAC();//config DAC  \n\n');
    fprintf(f,'\tconfigTimerInterrupt(%d);//config timer interrupt to be activate every observer sampling time  \n\n',object.getSamplingTime);
    
    fprintf(f,'\tz[0] = default_z;\n\n');
    
    fprintf(f,'\twhile(1)\n');
    fprintf(f,'\t{\n');
    fprintf(f,'\t\twhile(!controlGo);//wait for interrupt activation\n');
    fprintf(f,'\t\tcontrolGo = 0;\n');
    fprintf(f,'\t\treadInputADC(uyADC);\n');
    fprintf(f,'\t\tscaleInput(uyADC,uy);\n');
    fprintf(f,'\t\t/* uy = [u1,u2,...,y1,y2,...] */\n\n');
    
    fprintf(f,'\t\tif(count_buf_u != mu_VS)\n');
    fprintf(f,'\t\t\tcount_buf_u++;\n');
    fprintf(f,'\t\tif(count_buf_y != my_VS)\n');
    fprintf(f,'\t\t\tcount_buf_y++;\n');
    fprintf(f,'\t\tif(count_buf_z != mz_VS)\n');
    fprintf(f,'\t\t\tcount_buf_z++;\n');
    
    fprintf(f,'\t\tif(count_buf_u == mu_VS && count_buf_y == my_VS && count_buf_z == mz_VS)\n');
    fprintf(f,'\t\t\tfull_buf = 1;\n\n');
    
    if nu > 0
        fprintf(f,'\t\tfor(i=0;i<nu_VS;i++)\n');
        fprintf(f,'\t\t\tu[i] = uy[i];\n');
        fprintf(f,'\t\tfor(i=0;i<nu_VS;i++) {\n');
        fprintf(f,'\t\t\tfor(j=mu_VS-1;j>0;j--)\n');
        fprintf(f,'\t\t\t\tbuf_u[i][j] = buf_u[i][j-1];\n');
        fprintf(f,'\t\t\tbuf_u[i][0] = u[i];\n');
        fprintf(f,'\t\t}\n');
    end
    if ny > 0
        fprintf(f,'\t\tfor(i=0;i<ny_VS;i++)\n');
        fprintf(f,'\t\t\ty[i] = uy[nu_VS+i];\n');
        fprintf(f,'\t\tfor(i=0;i<ny_VS;i++) {\n');
        fprintf(f,'\t\t\tfor(j=my_VS-1;j>0;j--)\n');
        fprintf(f,'\t\t\t\tbuf_y[i][j] = buf_y[i][j-1];\n');
        fprintf(f,'\t\t\tbuf_y[i][0] = y[i];\n');
        fprintf(f,'\t\t}\n');
    end
    if mz > 0
        fprintf(f,'\t\tfor(i=0;i<nz_VS;i++) {\n');
        fprintf(f,'\t\t\tfor(j=mz_VS-1;j>0;j--)\n');
        fprintf(f,'\t\t\t\tbuf_z[i][j] = buf_z[i][j-1];\n');
        fprintf(f,'\t\t\tbuf_z[i][0] = z[i];\n');
        fprintf(f,'\t\t}\n');
    end
    
    fprintf(f,'\t\tif(full_buf == 1) {\n');
    fprintf(f,'\t\t\testimate(buf_u,buf_y,buf_z,z);\n');
    fprintf(f,'\t\t\tscaleOutput(zADC,z);\n');
    fprintf(f,'\t\t\tsetOutputDAC(zADC);\n');
    fprintf(f,'\t\t}\n');
    fprintf(f,'\t}\n');
    fprintf(f,'}\n\n');
    fprintf(f,'void scaleInput(int xADC[nu_VS+ny_VS], float xscale[nu_VS+ny_VS])\n{\n');
    
    k = 0;
    
    for i = 1:nu
        deltaReal =  domain.umax(i)-domain.umin(i);
        deltaADC = circuit_parameters.inputRange.max-circuit_parameters.inputRange.min;
        meanReal = (domain.umax(i)+domain.umin(i))/2;
        meanADC = (circuit_parameters.inputRange.max+circuit_parameters.inputRange.min)/2;
        fprintf(f,'\txscale[%d] = %.10e*xADC[%d]+%.10e;\n',k,deltaReal/deltaADC,k,(meanReal-deltaReal/deltaADC*meanADC));
        k = k+1;
    end
    
    for i = 1:ny
        deltaReal =  domain.ymax(i)-domain.ymin(i);
        deltaADC = circuit_parameters.inputRange.max-circuit_parameters.inputRange.min;
        meanReal = (domain.ymax(i)+domain.ymin(i))/2;
        meanADC = (circuit_parameters.inputRange.max+circuit_parameters.inputRange.min)/2;
        fprintf(f,'\txscale[%d] = %.10e*xADC[%d]+%.10e;\n',k,deltaReal/deltaADC,k,(meanReal-deltaReal/deltaADC*meanADC));
        k = k+1;
    end
    
    fprintf(f,'}\n\n');
    
    fprintf(f,'void scaleOutput(int uADC[nz_VS], float uscale[nz_VS])\n{\n');
    for i = 1:nz
        deltaReal =  domain.zmax(i)-domain.zmin(i);
        deltaADC = circuit_parameters.outputRange.max-circuit_parameters.outputRange.min;
        meanReal = (domain.zmax(i)+domain.zmin(i))/2;
        meanADC = (circuit_parameters.outputRange.max+circuit_parameters.outputRange.min)/2;
        
        fprintf(f,'\tuADC[%d] = %.10e*uscale[%d]+%.10e;\n',i-1,deltaADC/deltaReal,i-1,(meanADC-deltaADC/deltaReal*meanReal));
    end
    fprintf(f,'}\n');
    fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .H FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'virtualSensor.h'],'w');

fprintf(f,'#ifndef _VIRTUALSENSOR_H_\n');
fprintf(f,'#define _VIRTUALSENSOR_H_\n');
fprintf(f,'#define nu_VS %d\n',nu);
fprintf(f,'#define ny_VS %d\n',ny);
fprintf(f,'#define nz_VS %d\n',nz);
if object.isCurrent
    fprintf(f,'#define mu_VS %d\n',mu(1));
    fprintf(f,'#define my_VS %d\n',my(1));
else
    if mu(1) ~= 0
        fprintf(f,'#define mu_VS %d\n',mu(1)+1);
    else
        fprintf(f,'#define mu_VS 0\n');
    end
    if my(1) ~= 0
        fprintf(f,'#define my_VS %d\n',my(1)+1);
    else
        fprintf(f,'#define my_VS 0\n');
    end
end
fprintf(f,'#define mz_VS %d\n',mz(1));
fprintf(f,'#define default_z %f\n',circuit_parameters.initialCondition);

fprintf(f,'void estimate(float buf_u[nu_VS][mu_VS], float buf_y[ny_VS][my_VS], float buf_z[nz_VS][mz_VS], float *z);\n');
fprintf(f,'#endif\n');
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'virtualSensor.c'],'w');

fprintf(f,'#include "virtualSensor.h"\n');
for i = 1:npwas
    fprintf(f,'#include "pwas_fun_%d.h"\n',i);
end
fprintf(f,'\n');
fprintf(f,'void estimate(float buf_u[nu_VS][mu_VS], float buf_y[ny_VS][my_VS], float buf_z[nz_VS][mz_VS], float *z)\n{\n');
fprintf(f,'\tint i,j,k;\n');
for i = 1:npwas
    fprintf(f,'\tfloat x%d[%d];\n',i-1,fun(i).getDomainDimensions());
    fprintf(f,'\tfloat z%d[1];\n',i-1);
end
if ~reducedComplexity
    fprintf(f,'\tk = 0;\n');
    fprintf(f,'\tfor(i=0;i<nu_VS;i++){\n');
    fprintf(f,'\t\tfor(j=0;j<mu_VS;j++){\n');
    fprintf(f,'\t\t\tx0[k] = buf_u[i][j];\n');
    fprintf(f,'\t\t\tk++;\n');
    fprintf(f,'\t\t}\n');
    fprintf(f,'\t}\n');
    fprintf(f,'\tfor(i=0;i<ny_VS;i++){\n');
    fprintf(f,'\t\tfor(j=0;j<my_VS;j++){\n');
    fprintf(f,'\t\t\tx0[k] = buf_y[i][j];\n');
    fprintf(f,'\t\t\tk++;\n');
    fprintf(f,'\t\t}\n');
    fprintf(f,'\t}\n');
    fprintf(f,'\tfor(i=0;i<nz_VS;i++){\n');
    fprintf(f,'\t\tfor(j=0;j<mz_VS;j++){\n');
    fprintf(f,'\t\t\tx0[k] = buf_z[i][j];\n');
    fprintf(f,'\t\t\tk++;\n');
    fprintf(f,'\t\t}\n');
    fprintf(f,'\t}\n');
    fprintf(f,'\tcalculatePWAS_1(x0,z0);\n');
    fprintf(f,'\tz[0] = z0[0];\n');
else
    
    unames = cell(npwas,nu);
    ynames = cell(npwas,ny);
    znames = cell(npwas,1);
    
    if object.isCurrent
        for i = 1:nu
            for j = 1:mu(i)
                unames{j,i} = ['buf_u[',num2str(i-1),'][',num2str(j-1),']'];
            end
        end
        for i = 1:ny
            for j = 1:my(i)
                ynames{j,i} = ['buf_y[',num2str(i-1),'][',num2str(j-1),']'];
            end
        end
        for j = 2:mz+1
            znames{j,1} = ['buf_z[0][',num2str(j-2),']'];
        end
    else
        for i = 1:nu
            for j = 1:mu(i)
                unames{j,i} = ['buf_u[',num2str(i-1),'][',num2str(j),']'];
            end
        end
        for i = 1:ny
            for j = 1:my(i)
                ynames{j,i} = ['buf_y[',num2str(i-1),'][',num2str(j),']'];
            end
        end
        for j = 1:mz
            znames{j,1} = ['buf_z[0][',num2str(j-1),']'];
        end
        
    end
    
    for ii = 1:npwas
        
        k = 0;
        for i = 1:nu
            if ~isempty(unames{ii,1})
                fprintf(f,'\tx%d[%d] = %s;\n',ii-1,k,unames{ii,i});
                k = k+1;
            end
        end
        for i = 1:ny
            if ~isempty(ynames{ii,i})
                fprintf(f,'\tx%d[%d] = %s;\n',ii-1,k,ynames{ii,i});
                k = k+1;
            end
        end
        for j = 1:mz
            if ~isempty(znames{ii,1})
                fprintf(f,'\tx%d[%d] = %s;\n',ii-1,k,znames{ii,1});
                k = k+1;
            end
        end
        
    end     
    
    for i = 1:npwas
        fprintf(f,'\tcalculatePWAS_%d(x%d,z%d);\n',i,i-1,i-1);
    end
    fprintf(f,'\tz[0] = ');
    for i = 1:npwas-1
        fprintf(f,'z%d[0]+',i-1);
    end
    fprintf(f,'z%d[0];\n',npwas-1);
end
fprintf(f,'}\n');

fclose(f);

% Prepare the output sructure
if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    range.umin = domain.umin;
    range.umax = domain.umax;
    range.ymin = domain.ymin;
    range.ymax = domain.ymax;
    range.zmin = domain.zmin;
    range.zmax = domain.zmax;
    
    cir_range.umin = circuit_parameters.inputRange.min;
    cir_range.umax = circuit_parameters.inputRange.max;
    cir_range.ymin = circuit_parameters.inputRange.min;
    cir_range.ymax = circuit_parameters.inputRange.max;
    
    cir_range.zmin = circuit_parameters.outputRange.min;
    cir_range.zmax = circuit_parameters.outputRange.max;
    
    %     optO.latency = latency;
    %     optO.memory_size = memory_size;
    %     optO.multipliers = multipliers;
    %     optO.nbitmul = nbit_mul;
    optO.range = range;
    optO.cir_range = cir_range;
    varargout{1} = optO;
end





% TO DO
% report

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % memory_size = funReport.memory_size;
% % %
% % % if nargout > 1
% % %     error('Wrong number of outputs')
% % % elseif nargout == 1
% % %     optO.memory_size = memory_size;
% % %     optO.range = range;
% % %     varargout{1} = optO;
% % % end
% % %
% % %
% % % % if needed generate interface
% % % if circuit_parameters.generate_main
% % %     filename = strcat(folder,'C_report.log');
% % %
% % %     fout = fopen(filename, 'w');
% % %
% % %     nbit = circuit_parameters.inputResolution;
% % %     nbitout = circuit_parameters.outputResolution;
% % %
% % %     fprintf(fout,'-------------------------------------------------------------\n');
% % %     fprintf(fout,'|                Circuit information report                  |\n');
% % %     fprintf(fout,'-------------------------------------------------------------\n\n');
% % %     fprintf(fout,'INPUTS\n');
% % %     fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.inputResolution);
% % %     fprintf(fout,'\t - Range (model --> circuit):\n');
% % %     for i = 1:nx
% % %         cirmin = decimal2signed(circuit_parameters.inputRange.min(i),nbit,0);
% % %         cirmax = decimal2signed(circuit_parameters.inputRange.max(i),nbit,0);
% % %         fprintf(fout,'\t\t%s: [%f %f] --> x%d: [%s %s]\n',object.xnames{i},range.xmin(i),range.xmax(i),i,cirmin.bin,cirmax.bin);
% % %     end
% % %     for i = 1:np
% % %         cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx),nbit,0);
% % %         cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx),nbit,0);
% % %         fprintf(fout,'\t\t%s: [%f %f] --> p%d: [%s %s]\n',object.pnames{i},range.pmin(i),range.pmax(i),i,cirmin.bin,cirmax.bin);
% % %     end
% % %     for i = 1:nd
% % %         cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np),nbit,0);
% % %         cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np),nbit,0);
% % %         fprintf(fout,'\t\t%s: [%f %f] --> d%d: [%s %s]\n',object.dnames{i},range.dmin(i),range.dmax(i),i,cirmin.bin,cirmax.bin);
% % %     end
% % %     for i = 1:nxref
% % %         track = object.getTrackingVariable;
% % %         cirmin = decimal2signed(circuit_parameters.inputRange.min(i+nx+np+nd),nbit,0);
% % %         cirmax = decimal2signed(circuit_parameters.inputRange.max(i+nx+np+nd),nbit,0);
% % %         fprintf(fout,'\t\t%s_ref: [%f %f] --> xref%d: [%s %s]\n',object.xnames{track(i)},range.xrefmin(i),range.xrefmax(i),i,cirmin.bin,cirmax.bin);
% % %     end
% % %     fprintf(fout,'\nOUTPUTS\n');
% % %     fprintf(fout,'\t - Resolution: %d bits\n',circuit_parameters.outputResolution);
% % %     fprintf(fout,'\t - Range (model --> circuit):\n');
% % %     for i = 1:nu
% % %         ucirmin = decimal2signed(circuit_parameters.outputRange.min(i),nbitout,0);
% % %         ucirmax = decimal2signed(circuit_parameters.outputRange.max(i),nbitout,0);
% % %         fprintf(fout,'\t\t%s: [%f %f] --> u%d: [%s %s]\n',object.unames{i},range.umin(i),range.umax(i),i,ucirmin.bin,ucirmax.bin);
% % %     end
% % %
% % %     fprintf(fout,'MEMORY SIZE\n');
% % %     fprintf(fout,'\t - Number of cells = %d\n',memory_size);
% % %     fprintf(fout,'\t - Word size = %d bits\n',32);
% % %     fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
% % %     fclose(fout);
% % %     edit([folder ,'C_report.log'])



end