function writePackage(object,nbit,nbit_coeff,nbitout,sampling_latency,defaultOutput,folder)
% writePackage   Writes the VHDL code implementing block controller_package
%
% writePackage(OBJ,OPTS)
% This is a private method

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Mateo Lodi (matteo.lodi@edu.unige.it)
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
% Boston, MA  02111-1307  US

% Actual number of bits to represent input: if input representation is
% unsigned, one more bit is necessary because the multipliers are signed
% and therefore a '0' must be added as MSB of the input
nbit_eff = nbit;

% Idem for the output
nbit_coeff_eff = nbit_coeff;

% Number of bits for the multiplier
nbit_mul = max(nbit_eff,nbit_coeff_eff);

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;

% Open files
fin = fopen([getvhdlpath(),'explicitMPCctrl/controller_package.vhd'],'r');
fout = fopen([folder,'controller_package.vhd'],'w');

while 1
    
    % Read line from input file
    tline = fgetl(fin);
    
    % If there are no characters, break
    if ~ischar(tline)
        break
    end
    
    % If the MATLABGEN section starts, write the VHDL code ...
    if strcmp(tline,'--- BEGIN MATLABGEN ---')
        
        fprintf(fout,tline);
        fprintf(fout,'\n');
        fprintf(fout,'\n');
        
        fprintf(fout,'\t-- Number of domain dimensions\n');
        fprintf(fout,'\tconstant N_DIM_CTRL : integer := %d;\n',object.getNumberOfDimensions);
        fprintf(fout,'\t-- Number of codomain dimensions\n');
        fprintf(fout,'\tconstant N_FUN_CTRL : integer := %d;\n',object.getNumberOfInputs);
        fprintf(fout,'\t-- Number of bits to represent the inputs\n');
        fprintf(fout,'\tconstant N_BIT_CTRL : integer := %d;\n',nbit);
        fprintf(fout,'\t-- Number of bits to represent the outputs\n');
        fprintf(fout,'\tconstant N_BIT_OUT_CTRL : integer := %d;\n',nbitout);
        fprintf(fout,'\t-- Number of bits to represent the coefficients\n');
        fprintf(fout,'\tconstant N_BIT_COEFF_CTRL : integer := %d;\n',nbit_coeff);
        fprintf(fout,'\t-- Number of bits to represent memory address\n');
        fprintf(fout,'\tconstant N_BIT_MUL_CTRL : integer := %d;\n',nbit_mul);
        fprintf(fout,'\t-- Definition of new types\n');
        fprintf(fout,'\t\n');
        
        fprintf(fout,'\tconstant sampling_latency_CTRL : integer := %d;\n\n',fix(sampling_latency));
        
        fprintf(fout,'\tconstant nx_CTRL : integer := %d;\n',object.getNumberOfStates);
        fprintf(fout,'\tconstant npar_CTRL : integer := %d;\n',object.getNumberOfParameters);
        fprintf(fout,'\tconstant nd_CTRL : integer := %d;\n',object.getNumberOfUnmeasurableInputs);
        fprintf(fout,'\tconstant nu_CTRL : integer := %d;\n',object.getNumberOfInputs);
        fprintf(fout,'\tconstant nxref_CTRL : integer := %d;\n',numel(object.getTrackingVariable));
        
        if nx ~= 0
            fprintf(fout,'\ttype xType_CTRL is array(nx_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        end
        if nd ~= 0
            fprintf(fout,'\ttype dType_CTRL is array(nd_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        end
        if nu ~= 0
            fprintf(fout,'\ttype uType_CTRL is array(nu_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        end
        if np ~= 0
            fprintf(fout,'\ttype pType_CTRL is array(npar_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        end
        
        
        
        fprintf(fout,'\ttype x_matrix_CTRL is array(N_DIM_CTRL-1 downto 0) of signed(N_BIT_CTRL-1 downto 0);\n');
        
        fprintf(fout,'\ttype u_matrix_CTRL is array(N_FUN_CTRL-1 downto 0) of signed(N_BIT_OUT_CTRL-1 downto 0);\n');
        fprintf(fout,'\ttype mul_in_matrix_CTRL is array(N_DIM_CTRL-1 downto 0) of signed(N_BIT_MUL_CTRL-1 downto 0);\n');
        fprintf(fout,'\ttype mul_out_matrix_CTRL is array(N_DIM_CTRL-1 downto 0) of signed(2*N_BIT_MUL_CTRL-1 downto 0);\n');
         if nu == 1
            fprintf(fout,'\ttype u_matrix_CTRL_toy is array(1 downto 0) of signed(N_BIT_OUT_CTRL-1 downto 0);\n');
            
            fprintf(fout,'\tconstant defaultOutput : u_matrix_CTRL_toy := (');
            
            fprintf(fout,'"%s","%s");\n\n',defaultOutput{1},defaultOutput{1});
        else
            fprintf(fout,'\tconstant defaultOutput : u_matrix_CTRL := (');
            for i=nu:-1:2
                fprintf(fout,'"%s",',defaultOutput{i});
            end
            fprintf(fout,'"%s");\n\n',defaultOutput{1});
        end
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'controller_package.vhd']);