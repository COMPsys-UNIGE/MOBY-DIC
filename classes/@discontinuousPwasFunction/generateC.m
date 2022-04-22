%% generateC
% Generates C files for the circuit implementation of the
% PWAS function on microcontroller
% SYNTAX
%
% object = generateVHDL(object)
% object = generateVHDL(object,opts)
% object = generateVHDL(object,idx,opts)
%
%  opts is a structure with the following fields:
% * generate_main: if this field ha value true the function work as
%     generateC(OBJ), otherwise the function doesn't generate the main file
% * folder : folde
% idx is a vector containing the functions that must be generated (Default
%         all function are generated)
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

circuit_parameters = object.discontinuousPwasCset(circuit_parameters,idx);


nfun = object.getCodomainDimensions;
nDim = object.getDomainDimensions;
% Set default value for idx
if isempty(idx)
    idx = 1:nfun;
end

folder = circuit_parameters.folder;

nDyn = object.nDyn;

ndim = object.getDomainDimensions;
nY = object.getCodomainDimensions;

[Hd Kd] = object.getDomain;
P = Polyhedron(Hd,Kd);

if ~P.isBounded
    error(['PWAS function domain must be bounded in order'...
        ' to generate circuit file']);
end
Dp = P.outerApprox();
DVert = Dp.V;

% w = object.getWeights;

minV = min(DVert);
maxV = max(DVert);

range.xmin = minV;
range.xmax = maxV;

w = object.functions(1).pwasFunction.getWeights;
wmin = min(w);
wmax = max(w);

for i=2:nDyn
    w = object.functions(i).pwasFunction.getWeights;
    wmin = min([w;wmin]);
    wmax = max([w;wmax]);
end

range.umin = wmin(idx);
range.umax = wmax(idx);

np = object.getNumberOfPartitions;
P = object.getPartition;

[~,DD] = object.getDomain;


if ~exist(folder)
    mkdir(folder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE MAIN FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if circuit_parameters.generate_main
    f = fopen([folder,'main.c'],'w');
    fprintf(f,'#include <stdio.h>\n');
    fprintf(f,'#include "discontinuouspwas_fun.h"\n');
    fprintf(f,'\n\n');
    fprintf(f,'int findDynamic(float x[nDim]);\n');
    fprintf(f,'\n\n');
    fprintf(f,'int main()\n{\n');
    fprintf(f,'float x[nDim];//input vector\n');
    fprintf(f,'float f[nY];//vector containYg function value \n\n');
    fprintf(f,'int i;\n');
    fprintf(f,'/* Fill x for example from file or from std input as : */\n');
    fprintf(f,'/* x = [x1,x2,...,xn,p1,..,pq,d1,...,dl]               */\n\n');
    fprintf(f,'i = findDynamic(x);\n');
    fprintf(f,'calculatediscontinuousPWAS(x,f,i);\n');
    fprintf(f,'}\n\n');
    
    
    HH = [];
    KK = [];
    
    reg = cell(object.nDyn,1);
    ff = object.functions;
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
    fprintf(f,'int findDynamic(float x[nDim])\n{\n');
    fprintf(f,'const float Hedge[%d][nDim] = {{',numel(Kred));
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
    for i=1:numel(Kred)-1PWAS
        fprintf(f,'%.10e,\n',Kred(i));
    end
    fprintf(f,'%f};\n\n',Kred(end));
    fprintf(f,'int i,j,res;\n');
    fprintf(f,'float tmp;\n');
    fprintf(f,'res = 0;\n');
    fprintf(f,'for(i=0;i<%d;i++)\n{\n',numel(Kred));
    fprintf(f,'tmp = 0;\n');
    fprintf(f,'for(j=0;j<nDim;j++)\n');
    fprintf(f,'    tmp += Hedge[i][j]*x[j];\n');
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

fout = fopen([folder,'discontinuouspwas_fun.h'],'w');
fprintf(fout,'#ifndef _PWASD_FUN_H_\n');
fprintf(fout,'#define _PWASD_FUN_H_\n\n');

fprintf(fout,'#define nDim %d\n',ndim);
fprintf(fout,'#define nY %d\n',object.getCodomainDimensions);
fprintf(fout,'#define nDyn %d\n',nDyn);
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

fprintf(fout,'void calculatediscontinuousPWAS(float *x, float *u,int dynamic);\n\n');

fprintf(fout,'float uFunction(int *vertexAddr, float *mu, int uIndex,int dynamic, float *x); \n\n');


fprintf(fout,'#endif\n');


fclose(fout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE .C FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen([folder,'discontinuouspwas_fun.c'],'w');

fprintf(fout,'#include "discontinuouspwas_fun.h"\n');

fprintf(fout,'const int Npartition[nDim] = {');
for i=1:numel(np)
    fprintf(fout,'%d',np(i));
    if i ~= numel(np)
        fprintf(fout,',');
    else
        fprintf(fout,'};\n\n');
    end
end

% Write m and q used to transform domain
if object.isUniform
    fprintf(fout,'const float m[nDim] = {');
    for i=1:ndim
        fprintf(fout,'%.10e',np(i)/(DD(i)+DD(ndim+i)));
        if i ~= ndim
            fprintf(fout,',');
        end
    end
    fprintf(fout,'};\n');
    
    
    
    
    fprintf(fout,'const float q[nDim] = {');
    for i=1:ndim
        fprintf(fout,'%.10e',-np(i)/(DD(i)+DD(ndim+i))*(-DD(ndim+i)));
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
            fprintf(fout,'%.10e',1/pLen(j));
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
            fprintf(fout,'%.10e',PP(j));
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
fprintf(fout,'const float weigthVector[nDyn][nY][nWeight] = {');
for k=1:nDyn
    fprintf(fout,'{');
    w = object.functions(i).pwasFunction.getWeights;
    for i=idx
        fprintf(fout,'{');
        for j=1:size(w,1)
            fprintf(fout,'%.10e',w(j,i));
            if j ~= size(w,1)
                fprintf(fout,',');
            else
                fprintf(fout,'}');
            end
        end
        if i ~= idx(end)
            fprintf(fout,',\n');
        else
            fprintf(fout,'}\n');
        end
    end
    if k ~= nDyn
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
    fprintf(fout,'    z[k] = m[k] * x[k] + q[k];\n');
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
fprintf(fout,'void calculatediscontinuousPWAS(float *x, float *u,int dynamic)\n');
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
fprintf(fout,'		*(u+i) = uFunction(vertAddr,mu,i,dynamic,Z);\n');
fprintf(fout,'}\n\n');

fprintf(fout,'float uFunction(int *vertexAddr, float *mu, int uIndex,int dynamic, float *x) \n');
fprintf(fout,'{\n');
fprintf(fout,'	int i;\n');
fprintf(fout,'	float res = 0;\n');
fprintf(fout,'	for(i=0;i<nDim+1;i++)\n');
fprintf(fout,'	{\n');
fprintf(fout,'      if(mu[i] != 0)\n');
fprintf(fout,'          res += mu[i]*weigthVector[dynamic][uIndex][vertexAddr[i]];\n');
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
        fprintf(fout,'\t\t%s: [%f %f]\n',object.xnames{i},range.xmin(i),range.xmax(i));
    end
    fprintf(fout,'\nOUTPUTS\n');
    fprintf(fout,'\t - Range:\n');
    for i = idx
        fprintf(fout,'\t\t%s: [%f %f]\n',object.ynames{i},range.umin(i),range.umax(i));
    end
    
    fprintf(fout,'MEMORY SIZE\n');
    fprintf(fout,'\t - Number of cells = %d\n',memory_size);
    fprintf(fout,'\t - Word size = %d bits\n',32);
    fprintf(fout,'\t - Total occupation = %.3f bytes\n',32*memory_size/8);
    fclose(fout);
    edit([folder ,'C_report.log'])
    
    
    
end

end