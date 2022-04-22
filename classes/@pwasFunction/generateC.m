function varargout = generateC(object,varargin)
% generateC Generates C files for the circuit implementation of the PWAS
%           function on microcontroller
%
% generateC(OBJ)
% Generates the C files for the circuit implementation on microcontroller
% of all components of the possibly vector PWAS function with the default
% options. Three files are generated:
%  - pwas_fun.h: header file containing the definition of constants and the 
%                declaration of all functions;
%  - pwas_fun.c: C file containing the implementation of all functions;
%  - main.c: main file where it is shown how to call the function for the 
%            computation of the PWAS function.
% A report is also generated showing the main circuit features. 
% The data are represented as float, therefore a 32bit floating point data 
% representation is used.
%
% generateC(OBJ,IDX)
% Generates the C files for the circuit implementation on microcontroller
% of the components of the vector PWAS function indicated by IDX, with the
% default options. IDX is a vector containing the indices of the components
% to implement. 
%
% generateC(OBJ,OPTS)
% Generates the VHDL files for the circuit implementation on microcontroller
% of all components of the possibly vector PWAS function with custom options.
% OPTS is a structure with the following fields:
%  - generateMain: flag indicating if the main file must be generated. 
%                  If the PWAS function is part of a controller, the main 
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

optout = pwasCset(circuit_parameters);


nfun = object.getCodomainDimensions;
% Set default value for idx
if isempty(idx)
    idx = 1:nfun;
end

folder = optout.folder;

ndim = object.getDomainDimensions;

[Hd, Kd] = object.getDomain;
P = Polyhedron(Hd,Kd);

if ~P.isBounded
    error(['PWAS function domain must be bounded in order'...
        ' to generate circuit file']);
end
Dp = P.outerApprox();
DVert = Dp.V;

w = object.getWeights;

minV = min(DVert);
maxV = max(DVert);

range.xmin = minV;
range.xmax = maxV;

wmin = min(w);
wmax = max(w);

range.umin = wmin(idx);
range.umax = wmax(idx);

np = object.getNumberOfPartitions;
P = object.getPartition;

[~,DD] = object.getDomain;


if ~exist(folder,'dir')
    mkdir(folder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE MAIN FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if optout.generateMain
    f = fopen([folder,'main.c'],'w');
    fprintf(f,'#include <stdio.h>\n');
    fprintf(f,'#include "pwas_fun.h"\n');
    fprintf(f,'\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'\tfloat x[nDim]; //vector of function inputs\n');
    fprintf(f,'\tfloat f[nY];   //vector of function outputs\n\n');
    fprintf(f,'\t// Read vector x (e.g., from file or standard input)\n');
    fprintf(f,'\t// ...\n\n');
    fprintf(f,'\t// Call function for the computation of the PWAS fuction\n');
    fprintf(f,'\tcalculatePWAS(x,f);\n\n');
    fprintf(f,'\t// Write vector f (e.g., to file or standard output)\n');
    fprintf(f,'\t// ...\n\n');
    fprintf(f,'}');
    fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .H FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen([folder,'pwas_fun.h'],'w');
fprintf(fout,'#ifndef _PWAS_FUN_H_\n');
fprintf(fout,'#define _PWAS_FUN_H_\n\n');

fprintf(fout,'#define nDim %d\n',ndim);
fprintf(fout,'#define nY %d\n',object.getCodomainDimensions);
fprintf(fout,'#define nWeight %d\n\n',size(w,1));
if ~object.isUniform
    fprintf(fout,'#define totalPartiton %d\n\n',sum(np));
end

fprintf(fout,'void getMu(float *z_sorted, float *mu);\n\n');

fprintf(fout,'void Sorting(float *number, int n);\n\n');

fprintf(fout,'void Transform(float *x, float *z);\n\n');

if object.isUniform
    fprintf(fout,'void divideZ(float *Z, int *intZ, float *decZ);\n\n');
    fprintf(fout,'void getSimplexVertices(int *intZ, float *decZ, float *sortedDecZ, int *addr);\n\n');
else
    fprintf(fout,'void divideZ(float *Z,int *partIndex, float *decZ);\n\n');
    fprintf(fout,'void getSimplexVertices(int *partIndex, float *decZ, float *sortedDecZ, int *addr);\n\n');
    fprintf(fout,'void findPartition(float *Z,int *partIndex);\n\n');
end

fprintf(fout,'void calculatePWAS(float *x, float *u);\n\n');

fprintf(fout,'float uFunction(int *vertexAddr, float *mu, int uIndex, float *x); \n\n');


fprintf(fout,'#endif\n');


fclose(fout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen([folder,'pwas_fun.c'],'w');

fprintf(fout,'#include "pwas_fun.h"\n');

fprintf(fout,'const int Npartition[nDim] = {');
for i=1:numel(np)
    fprintf(fout,'%d',np(i));
    if i ~= numel(np)
        fprintf(fout,',');
    else
        fprintf(fout,'};\n\n');
    end
end

% Write mm and qq used to transform domain
if object.isUniform
    fprintf(fout,'const float mm[nDim] = {');
    for i=1:ndim
        fprintf(fout,'%f',np(i)/(DD(i)+DD(ndim+i)));
        if i ~= ndim
            fprintf(fout,',');
        end
    end
    fprintf(fout,'};\n');
    
    fprintf(fout,'const float qq[nDim] = {');
    for i=1:ndim
        fprintf(fout,'%f',-np(i)/(DD(i)+DD(ndim+i))*(-DD(ndim+i)));
        if i ~= ndim
            fprintf(fout,',');
        end
    end
    fprintf(fout,'};\n\n\n');
end
%Write lengh of each partition
if ~object.isUniform
    fprintf(fout,'//linear vector containing lenght of each partition\n');
    fprintf(fout,'const float one_over_partLenght[totalPartiton] = {');
    for i=1:ndim
        pLen = diff(P{i});
        for j=1:np(i)
            fprintf(fout,'%f',1/pLen(j));
            if j~= np(i)
                fprintf(fout,',');
            else
                if i~=ndim
                    fprintf(fout,',\n');
                end
            end
        end
    end
    fprintf(fout,'};\n');
    
    fprintf(fout,'//linear vector containing partiton left edge\n');
    fprintf(fout,'const float partLeftEdge[totalPartiton+nDim] = {');
    for i=1:ndim
        PP = P{i};
        for j=1:np(i)+1
            fprintf(fout,'%f',PP(j));
            if j~= np(i)+1
                fprintf(fout,',');
            else
                if i~=ndim
                    fprintf(fout,',\n');
                end
            end
        end
    end
    fprintf(fout,'};\n');
    
    
end

% Write weights
fprintf(fout,'const float weigthVector[nY][nWeight] = {');
for i=idx
    fprintf(fout,'{');
    for j=1:size(w,1)
        fprintf(fout,'%f',w(j,i));
        if j ~= size(w,1)
            fprintf(fout,',');
        else
            fprintf(fout,'}');
        end
    end
    if i ~= idx(end)
        fprintf(fout,',\n');
    else
        fprintf(fout,'};\n\n');
    end
end

fprintf(fout,'void getMu(float *z_sorted,float *mu)\n');
fprintf(fout,'{\n');
fprintf(fout,'	int i;\n');
fprintf(fout,'	mu[0] = 1 - z_sorted[0];\n\n');

fprintf(fout,'	for (i = 0; i < nDim; i++)\n');
fprintf(fout,'	{\n');
fprintf(fout,'			mu[i + 1] = z_sorted[i] - z_sorted[i + 1];\n');
fprintf(fout,'	}\n');
fprintf(fout,'	mu[nDim] = z_sorted[nDim-1];\n');
fprintf(fout,'}\n\n');

fprintf(fout,'void Sorting(float *number, int n)\n');
fprintf(fout,'{\n');
fprintf(fout,'	int i,j;\n');
fprintf(fout,'	float tmp;\n\n');
fprintf(fout,'	for (i = 0; i < n; ++i)\n');
fprintf(fout,'	{\n');
fprintf(fout,'		for (j = i + 1; j < n; ++j)\n');
fprintf(fout,'		{\n');
fprintf(fout,'			if (number[i] < number[j])\n');
fprintf(fout,'			{\n');
fprintf(fout,'				tmp = number[i];\n');
fprintf(fout,'				number[i] = number[j];\n');
fprintf(fout,'				number[j] = tmp;\n');
fprintf(fout,'			}\n');
fprintf(fout,'		}\n');
fprintf(fout,'	}\n');
fprintf(fout,'}\n\n');

if object.isUniform
    fprintf(fout,'void Transform(float *x, float *z)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int k;\n');
    fprintf(fout,'	int i;\n');
    fprintf(fout,'	for (k = 0; k < nDim; k++)\n');
    fprintf(fout,'	{\n');
    fprintf(fout,'    z[k] = mm[k] * x[k] + qq[k];\n');
    fprintf(fout,'	}\n');
    fprintf(fout,'}\n\n');
end

if object.isUniform
    fprintf(fout,'void divideZ(float *Z, int *intZ, float *decZ)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int i;\n\n');
    fprintf(fout,'	for (i = 0; i < nDim; i++)\n');
    fprintf(fout,'	{\n');
    fprintf(fout,'		intZ[i] = (int)(Z[i]);\n');
    fprintf(fout,'		decZ[i] = Z[i] - intZ[i];\n');
    fprintf(fout,'	}\n');
    fprintf(fout,'}\n\n');
else
    fprintf(fout,'void divideZ(float *Z,int *partIndex, float *decZ)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int i;\n\n');
    fprintf(fout,'  int pastVectLen = 0;\n');
    fprintf(fout,'	for (i = 0; i < nDim; i++)\n');
    fprintf(fout,'	{\n');
    fprintf(fout,'		decZ[i] = (Z[i] - partLeftEdge[partIndex[i]+i+pastVectLen])*one_over_partLenght[pastVectLen+partIndex[i]];\n');
    fprintf(fout,'  	pastVectLen += Npartition[i];\n');
    fprintf(fout,'	}\n');
    fprintf(fout,'}\n\n');
end


if object.isUniform
    fprintf(fout,'void getSimplexVertices(int *intZ, float *decZ, float *sortedDecZ, int *addr)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int i,j;\n');
    fprintf(fout,'	int toAdd;\n');
    fprintf(fout,'	float diff;\n');
    fprintf(fout,'	addr[0] = 0;\n');
    fprintf(fout,'	for(i=nDim-1;i>0;i--)\n');
    fprintf(fout,'		addr[0] = (addr[0]+intZ[i])*(Npartition[i-1]+1);\n');
    fprintf(fout,'	addr[0] += intZ[0];\n\n');
else
    fprintf(fout,'void getSimplexVertices(int *partIndex, float *decZ, float *sortedDecZ, int *addr)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int i,j;\n');
    fprintf(fout,'	int toAdd;\n');
    fprintf(fout,'	float diff;\n');
    fprintf(fout,'	addr[0] = 0;\n');
    fprintf(fout,'	for(i=nDim-1;i>0;i--)\n');
    fprintf(fout,'		addr[0] = (addr[0]+partIndex[i])*(Npartition[i-1]+1);\n');
    fprintf(fout,'	addr[0] += partIndex[0];\n\n');
end
fprintf(fout,'	for(i=1;i<nDim+1;i++)\n');
fprintf(fout,'	{\n');
fprintf(fout,'		addr[i] = 0;\n');
fprintf(fout,'		for(j=nDim-1;j>0;j--)\n');
fprintf(fout,'		{\n');
fprintf(fout,'			//1(decZ[j]-sortedDecZ[i]), where 1 is the step function\n');
fprintf(fout,'			diff = decZ[j]-sortedDecZ[i-1]; \n');
fprintf(fout,'			if (diff >= 0)\n');
fprintf(fout,'				toAdd = 1;\n');
fprintf(fout,'			else \n');
fprintf(fout,'				toAdd = 0;\n');
if object.isUniform
    fprintf(fout,'		addr[i] = (addr[i]+intZ[j]+toAdd)*(Npartition[j-1]+1);\n');
else
    fprintf(fout,'		addr[i] = (addr[i]+partIndex[j]+toAdd)*(Npartition[j-1]+1);\n');
end
fprintf(fout,'		}\n');
fprintf(fout,'	diff = decZ[0]-sortedDecZ[i-1]; \n');
fprintf(fout,'		if (diff >= 0)\n');
fprintf(fout,'			toAdd = 1;\n');
fprintf(fout,'		else \n');
fprintf(fout,'			toAdd = 0;\n');
if object.isUniform
    fprintf(fout,'      	addr[i] += intZ[0]+toAdd;\n');
else
    fprintf(fout,'      	addr[i] += partIndex[0]+toAdd;\n');
end
fprintf(fout,'	}\n');
fprintf(fout,'}\n\n');
if ~object.isUniform
    fprintf(fout,'void findPartition(float *Z,int *partIndex)\n');
    fprintf(fout,'{\n');
    fprintf(fout,'	int i,j;\n');
    fprintf(fout,'  int pastVectLen = 0;\n');
    fprintf(fout,'	for(i=0;i<nDim;i++)\n');
    fprintf(fout,'  {\n');
    fprintf(fout,'      for(j=0;j<Npartition[i]+1;j++)\n');
    fprintf(fout,'      {\n');
    fprintf(fout,'          if (Z[i] < partLeftEdge[pastVectLen+j])\n');
    fprintf(fout,'          	break;\n');
    fprintf(fout,'      }\n');
    fprintf(fout,'   partIndex[i] = j-1;\n');
    fprintf(fout,'   pastVectLen += Npartition[i]+1;\n');
    fprintf(fout,'   }\n');
    fprintf(fout,'}\n\n');
end
fprintf(fout,'void calculatePWAS(float *x, float *u)\n');
fprintf(fout,'{\n');
fprintf(fout,'	int i;\n\n');

if object.isUniform
    fprintf(fout,'	float Z[nDim];\n');
    fprintf(fout,'	int Zint[nDim];\n');
else
    fprintf(fout,'	float *Z;\n');
    fprintf(fout,'	int partIndex[nDim];\n');
end

fprintf(fout,'	float Zdec[nDim];\n');
fprintf(fout,'	float ZdecSorted[nDim];\n\n');

fprintf(fout,'	float mu[nDim+1];\n');
fprintf(fout,'	int vertAddr[nDim+1];\n\n');
if object.isUniform
    fprintf(fout,'	Transform(x, Z);\n');
    fprintf(fout,'	divideZ(Z, Zint, Zdec);\n');
else
    fprintf(fout,'	Z = x;\n');
    fprintf(fout,'  findPartition(Z,partIndex);\n');
    fprintf(fout,'	divideZ(Z, partIndex, Zdec);\n');
end
fprintf(fout,'	for(i=0;i<nDim;i++)\n');
fprintf(fout,'	   ZdecSorted[i] = Zdec[i];\n');
fprintf(fout,'	Sorting(ZdecSorted,nDim);\n');
fprintf(fout,'	getMu(ZdecSorted,mu);\n');
if object.isUniform
    fprintf(fout,'	getSimplexVertices(Zint, Zdec, ZdecSorted, vertAddr);\n\n');
else
    fprintf(fout,'	getSimplexVertices(partIndex, Zdec, ZdecSorted, vertAddr);\n\n');
end

fprintf(fout,'	for(i=0;i<nY;i++)\n');
fprintf(fout,'		*(u+i) = uFunction(vertAddr,mu,i,Z);\n');
fprintf(fout,'}\n\n');

fprintf(fout,'float uFunction(int *vertexAddr, float *mu, int uIndex, float *x) \n');
fprintf(fout,'{\n');
fprintf(fout,'	int i;\n');
fprintf(fout,'	float res = 0;\n');
fprintf(fout,'	for(i=0;i<nDim+1;i++)\n');
fprintf(fout,'	{\n');
fprintf(fout,'      if(mu[i] != 0)\n');
fprintf(fout,'          res += mu[i]*weigthVector[uIndex][vertexAddr[i]];\n');
fprintf(fout,'	}\n');
fprintf(fout,'	return res;\n');
fprintf(fout,'}\n');

fclose(fout);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memory_size = size(w,1)*numel(idx);

if nargout > 1
    error('Wrong number of outputs')
elseif nargout == 1
    optO.memory_size = memory_size;
    optO.range = range;
    varargout{1} = optO;
end

filename = strcat(folder,'C_report.log');

fout = fopen(filename, 'w');

fprintf(fout,'                        CIRCUIT INFORMATION REPORT\n');
fprintf(fout,'___________________________________________________________________________\n\n');
fprintf(fout,'Target device: microcontroller\n');
fprintf(fout,'___________________________________________________________________________\n\n');
fprintf(fout,'Inputs:\n');
fprintf(fout,'\t - Range (model):\n');
for i = 1:ndim
    fprintf(fout,'\t\t%s: [%f %f]\n',object.xnames{i},range.xmin(i),range.xmax(i));
end
fprintf(fout,'\nOutputs:\n');
fprintf(fout,'\t - Range (model):\n');
for i = 1:nfun
    fprintf(fout,'\t\t%s: [%f %f]\n',object.ynames{idx(i)},range.umin(i),range.umax(i));
end
fprintf(fout,'\nCoefficients:\n');
fprintf(fout,'\t - Representation: floating point signed\n');
fprintf(fout,'___________________________________________________________________________\n\n');
fprintf(fout,'Memory size:\n');
fprintf(fout,'\t - Number of cells = %d\n',memory_size);
fprintf(fout,'\t - Word size = 32 bits\n');
fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
fclose(fout);

if nargout == 0
    edit([folder ,'C_report.log'])
end

end

