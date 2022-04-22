function varargout =  generateC(object,varargin)
% generateC Generates C files for the circuit implementation of the PWAG
%           function on microcontroller
%
% generateC(OBJ)
% Generates the C files for the circuit implementation on microcontroller
% of all components of the possibly vector PWAG function with the default
% options. Three files are generated:
%  - pwag_fun.h: header file containing the definition of constants and the 
%                declaration of all functions;
%  - pwag_fun.c: C file containing the implementation of all functions;
%  - main.c: main file where it is shown how to call the function for the 
%            computation of the PWAG function.
% A report is also generated showing the main circuit features. 
% The data are represented as float, therefore a 32bit floating point data 
% representation is used.
%
% generateC(OBJ,IDX)
% Generates the C files for the circuit implementation on microcontroller
% of the components of the vector PWAG function indicated by IDX, with the
% default options. IDX is a vector containing the indices of the components
% to implement. 
%
% generateC(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation on microcontroller
% of all components of the possibly vector PWAG function with custom options.
% OPTS is a structure with the following fields:
%  - generateMain: flag indicating if the main file must be generated. 
%                  If the PWAG function is part of a controller, the main 
%                  file must not be generated. Default value: true.
%  - folder: folder where C files are created. Default value: progressive folder.
%
% PERF = generateC(OBJ,...)
% If an output argument is provided, the report is not automatically opened
% and a structure PERF containing some information about the circuit 
% implementation is returned. PERF is a structure with the following fields:
%  - memory_size: number of elements in the circuit memory
%  - range: structure with the following fields:
%            - xmin, xmax: ranges of the function inputs
%            - umin, umax: ranges of the function outputs

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


nfun = object.getCodomainDimensions;
nDim = object.getDomainDimensions;
% Set default value for idx
if isempty(idx)
    idx = 1:nfun;
end


circuit_parameters = object.pwagCset(circuit_parameters,idx);

folder = circuit_parameters.folder;

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

constr = object.getDomain;

umin = repmat(Inf,nfun,1);
umax = repmat(-Inf,nfun,1);

for i=1:object.getNumberOfRegions
    r = object.getRegions(i);
    Ri = Polyhedron(r.H,r.K);
    vertici = Ri.V;
    for i=1:nfun
        umin(i) = min([umin(i),r.F(i,:)*vertici'+r.G(i)]);
        umax(i) = max([umax(i),r.F(i,:)*vertici'+r.G(i)]);
    end
end

range.umin = umin(idx)';
range.umax = umax(idx)';

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

nNodes = numel(stateMem);

if ~exist(folder)
    mkdir(folder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE MAIN FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if circuit_parameters.generate_main
f = fopen([folder,'main.c'],'w');
fprintf(f,'#include <stdio.h>\n');
fprintf(f,'#include "pwag_fun.h"\n');
fprintf(f,'\n\n');
fprintf(f,'int main()\n{\n');
fprintf(f,'float x[nDim];//input vector\n');
fprintf(f,'float f[nY];//vector containYg function value \n\n');
fprintf(f,'/* Fill x for example from file or from std input as : */\n');
fprintf(f,'/* x = [x1,x2,...,xn,p1,..,pq,d1,...,dl]               */\n\n');
fprintf(f,'calculatePWAG(x,f);\n');
fprintf(f,'}');
fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .H FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'pwag_fun.h'],'w');

fprintf(f,'#ifndef _PWAG_FUN_H_\n');
fprintf(f,'#define _PWAG_FUN_H_\n');

fprintf(f,'#define nDim %d\n',nDim);
fprintf(f,'#define nY %d\n',nu);
fprintf(f,'#define nNodes %d\n',numel(stateMem));
fprintf(f,'#define rowEdgeMatrix %d\n',numel(Kred));
fprintf(f,'#define rowFunctionMatrix %d\n',numel(Gred));
fprintf(f,'\n');


fprintf(f,'typedef struct node {\n');
fprintf(f,'\tchar leaf;\n');
fprintf(f,'\tint index_in_matrix[nY]; //indice di H e K se il nodo non è un nodo foglio o F e G altrimenti\n');
fprintf(f,'\tint nodeLeftIndex;	//Left Node\n');
fprintf(f,'\tint nodeRightIndex;	//Right Node\n');
fprintf(f,'} Node;\n');
fprintf(f,'\n\n');

fprintf(f,'void calculatePWAG(float *x, float *u);\n');
fprintf(f,'float uFunction(int FandG_Index, float *x) ;\n');
fprintf(f,'void searchRegionFromPoint(float *x, int ind[nY]);\n\n');

fprintf(f,'#endif\n');
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen([folder,'pwag_fun.c'],'w');

fprintf(f,'#include "pwag_fun.h"\n\n');

fprintf(f,'const  Node nodeArray[nNodes] = {');
for i=1:nNodes
    current_node = stateMem(i);
    isLeaf = current_node.leaf;
    if isLeaf == 1
        ind = current_node.addrFG-numel(Kred);
        leftChildren = -1;
        rightChildren = -1;
    else
        ind = [current_node.addrHK, repmat(-1,1,nu-1)];
        leftChildren = current_node.next(1)-1;
        rightChildren = current_node.next(2)-1;
    end
    
    fprintf(f,'{%d,{',isLeaf);
    for j=1:nu-1
        fprintf(f,'%d,',ind(j)-1);
    end
    fprintf(f,'%d},',ind(end)-1);
    fprintf(f,'%d,%d}',leftChildren,rightChildren);
    if i ~= nNodes
        fprintf(f,',\n\t\t');
    else
        fprintf(f,'};\n');
    end
end
fprintf(f,'\n');
% Write Matrices
% H
fprintf(f,'const  float H[rowEdgeMatrix][nDim] = {');
for j=1:numel(Kred)
    fprintf(f,'{');
    for k=1:nDim
        fprintf(f,'%d',Hred(j,k));
        if k ~= nDim
            fprintf(f,',');
        else
            fprintf(f,'}');
        end
    end
    if j ~= numel(Kred)
        fprintf(f,',\n\t');
    else
        fprintf(f,'};\n');
    end
end
fprintf(f,'\n');

%K
fprintf(f,'const  float K[rowEdgeMatrix] = {');
for j=1:numel(Kred)
    fprintf(f,'%d',Kred(j,1));
    if j ~= numel(Kred)
        fprintf(f,',\n\t\t');
    else
        fprintf(f,'};\n');
    end
end
fprintf(f,'\n');


%.10e
fprintf(f,'const  float F[rowFunctionMatrix][nDim] = {');
for i=1:numel(Gred)
    fprintf(f,'{');
    for j=1:nDim
            fprintf(f,'%d',Fred(i,j));
            if j ~= nDim
                fprintf(f,',');
            else
                fprintf(f,'}');
            end
    end
    if i ~= numel(Gred)
        fprintf(f,',\n\t');
    else
        fprintf(f,'};\n');
    end
end
fprintf(f,'\n');

%G
fprintf(f,'const  float G[rowFunctionMatrix] = {');
for i=1:numel(Gred)-1
    fprintf(f,'%d,\n',Gred(i));
end
 fprintf(f,'%d};',Gred(end));
fprintf(f,'\n\n');

fprintf(f,'void calculatePWAG(float *x, float *u)\n');
fprintf(f,'{\n');
fprintf(f,'\tint i;\n');
fprintf(f,'\tint regIndex[nY];\nsearchRegionFromPoint(x,regIndex);\n\n');
	
fprintf(f,'\tfor(i=0;i<nY;i++)\n');
fprintf(f,'\t\t*(u+i) = uFunction(regIndex[i],x);\n');
fprintf(f,'}\n');


fprintf(f,'float uFunction(int FandG_Index, float *x) \n');
fprintf(f,'{\n');
fprintf(f,'\tint j;\n');
fprintf(f,'\tfloat sum = 0;\n');
fprintf(f,'\tfor (j = 0; j < nDim; j++)\n');
fprintf(f,'\t\tsum += F[FandG_Index][j] * x[j];\n');
fprintf(f,'\tsum += G[FandG_Index];\n');
fprintf(f,'\treturn sum;\n');
fprintf(f,'}\n\n');

fprintf(f,'void searchRegionFromPoint(float *x,int ind[nY]) \n');
fprintf(f,'{\n');
fprintf(f,'//Start from root node\n');
fprintf(f,'int i;\n');
fprintf(f,'Node n = nodeArray[0];\n');
fprintf(f,'while (!n.leaf) \n');
fprintf(f,'\t{\n');
fprintf(f,'\tint j;\n');
fprintf(f,'\tint ind = n.index_in_matrix[0];\n');
fprintf(f,'\tfloat sum = 0;\n');
fprintf(f,'\tfor (j = 0; j < nDim; j++)\n');
fprintf(f,'\t\tsum += H[ind][j] * x[j];\n');
fprintf(f,'\tif (sum > K[ind])\n');
fprintf(f,'\tn = nodeArray[n.nodeLeftIndex];\n');
fprintf(f,'\telse\n');
fprintf(f,'\t\tn = nodeArray[n.nodeRightIndex];\n');
fprintf(f,'\t}\n');
fprintf(f,'for(i=0;i<nY;i++)\n');
fprintf(f,'\tind[i] = n.index_in_matrix[i];\n');
fprintf(f,'}\n\n');

fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memory_size = numel([Hred Kred;Fred Gred]);

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
    
    
    fprintf(fout,'-------------------------------------------------------------\n');
    fprintf(fout,'|                Circuit information report                  |\n');
    fprintf(fout,'-------------------------------------------------------------\n\n');
    fprintf(fout,'INPUTS\n');
    fprintf(fout,'\t - Range:\n');
    for i = 1:nDim
        fprintf(fout,'\t\t%s: [%.10e %.10e]\n',object.xnames{i},range.xmin(i),range.xmax(i));
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Range:\n');
    for i = idx
        fprintf(fout,'\t\t%s: [%.10e %.10e]\n',object.ynames{i},range.umin(i),range.umax(i));
    end
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])
    
    
    
end

end