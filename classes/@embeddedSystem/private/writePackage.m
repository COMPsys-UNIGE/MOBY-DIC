function writePackage(object,nbit,nbit_coeff,nbitout,nMul,nbitMul,nMul_ctrl,nMul_obs,obs_latency,ctrl_latency,nDyn,folder)
% writePackage
% Writes the VHDL code implementing package pwasFunctionPackage
%
% SYNTAX
% This is a private method
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Alberto Oliveri (alberto.oliveri@unige.it)
% * Tomaso Poggi (tpoggi@essbilbao.org)
%
% Copyright is with:
%
% * Copyright (C) 2010-2011 University of Genoa, Italy.

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


% Actual number of bits to represent input: if input representation is
% unsigned, one more bit is necessary because the multipliers are signed
% and therefore a '0' must be added as MSB of the input
nbit_eff = nbit;

% Idem for the output
nbit_coeff_eff = nbit_coeff;


% Number of bits for the multiplier
nbit_mul = nbitMul;

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;
nxref = numel(object.getController().getTrackingVariable());

% Open files
fin = fopen([getvhdlpath(),'embeddedSystem/embeddedSystemPackage.vhd'],'r');
fout = fopen([folder,'embeddedSystemPackage.vhd'],'w');

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
        
        %         fprintf(fout,'\t-- Number of domain dimensions\n');
        %         fprintf(fout,'\tconstant N_DIM_ES : integer := %d;\n',object.getNumberOfDimensions);
        %         fprintf(fout,'\t-- Number of codomain dimensions\n');
        %         fprintf(fout,'\tconstant N_FUN_ES : integer := %d;\n',object.getNumberOfInputs);
        %         fprintf(fout,'\t-- Number of bits to represent the inputs\n');
        fprintf(fout,'\tconstant N_BIT_ES : integer := %d;\n',nbit);
        fprintf(fout,'\t-- Number of bits to represent the outputs\n');
        fprintf(fout,'\tconstant N_BIT_OUT_ES : integer := %d;\n',nbitout);
        fprintf(fout,'\t-- Number of bits to represent the coefficients\n');
        fprintf(fout,'\tconstant N_BIT_COEFF_ES : integer := %d;\n',nbit_coeff);
        fprintf(fout,'\t-- Number of bits of the multiplier\n');
        fprintf(fout,'\tconstant N_BIT_MUL_ES : integer := %d;\n',nbit_mul);
        fprintf(fout,'\t-- Number of multiplier\n');
        fprintf(fout,'\tconstant N_MUL_ES : integer := %d;\n',nMul);
        fprintf(fout,'\t-- Number of multiplier used by controller\n');
        fprintf(fout,'\tconstant N_MUL_CTRL_ES : integer := %d;\n',nMul_ctrl);
        fprintf(fout,'\t-- Number of multiplier used by observer\n');
        fprintf(fout,'\tconstant N_MUL_OBS_ES : integer := %d;\n',nMul_obs);
        fprintf(fout,'\t-- Number of inputs \n');
        fprintf(fout,'\tconstant N_INPUT_ES : integer := %d;\n',np+ny+nxref);
        fprintf(fout,'\t-- Number of outputr\n');
        fprintf(fout,'\tconstant N_OUTPUT_ES : integer := %d;\n\n',nu);
        
        fprintf(fout,'\tconstant nx_ES : integer := %d;\n',nx);
        fprintf(fout,'\tconstant np_ES : integer := %d;\n',np);
        fprintf(fout,'\tconstant nd_ES : integer := %d;\n',nd);
        fprintf(fout,'\tconstant nu_ES : integer := %d;\n',nu);
        fprintf(fout,'\tconstant ny_ES : integer := %d;\n',ny);
        fprintf(fout,'\tconstant nxref_ES : integer := %d;\n\n',nxref);
        if ~isempty(nDyn)
                fprintf(fout,'\tconstant nDyn_ES : integer := %d;\n\n',nDyn);
        end
        
        fprintf(fout,'\tconstant observer_latency_ES : integer := %d;\n',fix(obs_latency));
        fprintf(fout,'\tconstant controller_latency_ES : integer := %d;\n',fix(ctrl_latency));
        fprintf(fout,'\tconstant controller_mul_observer_latency_ES : integer := %d;\n\n',fix(ctrl_latency/obs_latency));
        fprintf(fout,'\t-- Definition of new types\n');
        fprintf(fout,'\t\n');
        
        fprintf(fout,'\t-- [p;x_ref;y] \n');
        fprintf(fout,'\ttype in_matrix_ES is array(N_INPUT_ES-1 downto 0) of std_logic_vector(N_BIT_ES-1 downto 0);\n');
        fprintf(fout,'\t-- [u] \n');
        fprintf(fout,'\ttype out_matrix_ES is array(N_OUTPUT_ES-1 downto 0) of std_logic_vector(N_BIT_OUT_ES-1 downto 0);\n');
        
        %         fprintf(fout,'\tconstant nx_CTRL : integer := %d;\n',object.getNumberOfStates);
        %         fprintf(fout,'\tconstant npar_CTRL : integer := %d;\n',object.getNumberOfParameters);
        %         fprintf(fout,'\tconstant nd_CTRL : integer := %d;\n',object.getNumberOfUnmeasurableInputs);
        %         fprintf(fout,'\tconstant nu_CTRL : integer := %d;\n',object.getNumberOfInputs);
        %         fprintf(fout,'\tconstant nxref_CTRL : integer := %d;\n',numel(object.getTrackingVariable));
        %
        %         if nx ~= 0
        %             fprintf(fout,'\ttype xType_CTRL is array(nx_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        %         end
        %         if nd ~= 0
        %             fprintf(fout,'\ttype dType_CTRL is array(nd_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        %         end
        %         if nu ~= 0
        %             fprintf(fout,'\ttype uType_CTRL is array(nu_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        %         end
        %         if np ~= 0
        %             fprintf(fout,'\ttype pType_CTRL is array(npar_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);\n');
        %         end
        %
        %
        %         fprintf(fout,'\ttype x_matrix_CTRL is array(N_DIM_CTRL-1 downto 0) of signed(N_BIT_CTRL-1 downto 0);\n');
        %
        %         fprintf(fout,'\ttype u_matrix_CTRL is array(N_FUN_CTRL-1 downto 0) of signed(N_BIT_OUT_CTRL-1 downto 0);\n');
        fprintf(fout,'\ttype mul_in_matrix_ES is array(N_MUL_ES-1 downto 0) of signed(N_BIT_MUL_ES-1 downto 0);\n');
        fprintf(fout,'\ttype mul_out_matrix_ES is array(N_MUL_ES-1 downto 0) of signed(2*N_BIT_MUL_ES-1 downto 0);\n');
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'embeddedSystemPackage.vhd']);