function writePackage(object,idx,opts,bits,scaleInput,Pbin)
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

% Extract input arguments
folder = opts.folder;
nbit = bits.nbit;
nint = bits.nint;
nbit_coeff = bits.nbit_coeff;
nbitout = bits.nbitout;
inputRepresentation = opts.inputRepresentation;
outputRepresentation = opts.outputRepresentation;
outputRange = opts.outputRange;

% Number of domain and codomain dimensions
ndim = object.getDomainDimensions();
nfun = numel(idx);

nDyn = object.nDyn;

% Number of subdivisions per dimension
np = object.getNumberOfPartitions();

uniform = object.isUniform();

% Actual number of bits to represent input: if input representation is
% unsigned, one more bit is necessary because the multipliers are signed
% and therefore a '0' must be added as MSB of the input
if strcmpi(inputRepresentation,'signed')
    nbit_eff = nbit;
else
    nbit_eff = nbit+1;
end
% Idem for the output
if strcmpi(outputRepresentation,'signed')
    nbit_coeff_eff = nbit_coeff;
else
    nbit_coeff_eff = nbit_coeff+1;
end

% Number of bits for the multiplier

nbit_mul = max(nbit+1,nbit_coeff_eff);


% Scaling factors in binary representation
if uniform
    for i = 1:ndim
        b(i) = decimal2signed(scaleInput.b(i),nbit_eff,0);
        A(i) = decimal2signed(scaleInput.A(i),nbit_coeff);
        App(i) = A(i).pointposition;
        npbin(i) = decimal2signed(np(i)-2^(-App(i)),2*nbit_mul,App(i));
    end
else  
    P = object.getPartition();
    for i = 1:ndim
        App_tmp = zeros(numel(P{i})-1,1);
        for j = 1:numel(P{i})-1
            A_tmp(j) = decimal2signed(scaleInput.A{i}(j),nbit_coeff);
            App_tmp(j) = A_tmp(j).pointposition;
        end
        App(i) = min(App_tmp);
        for j = 1:numel(P{i})-1
            A{i}(j) = decimal2signed(scaleInput.A{i}(j),nbit_coeff,App(i));
            b{i}(j) = decimal2signed(scaleInput.b{i}(j),nbit_eff,0);
        end
        npbin(i) = decimal2signed(1-2^(-App(i)),2*nbit_mul,App(i));
    end
end

% Find the number of decimal digits for the output
for i = 1:nfun
    fmin = outputRange.min(i);
    fmax = outputRange.max(i);
    if strcmpi(outputRepresentation,'signed')
        fminb = decimal2signed(fmin,nbit_coeff);
        fmaxb = decimal2signed(fmax,nbit_coeff);
    else
        fminb = decimal2unsigned(fmin,nbit_coeff);
        fmaxb = decimal2unsigned(fmax,nbit_coeff);
    end
    fpp(i) = min(fminb.pointposition,fmaxb.pointposition);
end

% Open files
fin = fopen([getvhdlpath(),'discontinuousPwasFunction/discontinuousPwasFunctionPackage.vhd'],'r');
fout = fopen([folder,'discontinuousPwasFunctionPackage.vhd'],'w');

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
        fprintf(fout,'\tconstant N_DIM_PWAS : integer := %d;\n',ndim);
        fprintf(fout,'\t-- Number of dynamic\n');
        fprintf(fout,'\tconstant N_DYN_PWAS : integer := %d;\n',nDyn);
        fprintf(fout,'\t-- Number of codomain dimensions\n');
        fprintf(fout,'\tconstant N_FUN_PWAS : integer := %d;\n',nfun);       
        fprintf(fout,'\t-- Number of bits to represent the inputs\n');
        fprintf(fout,'\tconstant N_BIT_PWAS : integer := %d;\n',nbit);
        fprintf(fout,'\t-- Number of bits to represent the outputs\n');
        fprintf(fout,'\tconstant N_BIT_OUT_PWAS : integer := %d;\n',nbitout);
        fprintf(fout,'\t-- Number of bits to represent the coefficients\n');
        fprintf(fout,'\tconstant N_BIT_COEFF_PWAS : integer := %d;\n',nbit_coeff);
        fprintf(fout,'\t-- Maximum number of bits used to represent the integer part of the scaled inputs\n');
        fprintf(fout,'\tconstant N_INT_MAX_PWAS : integer := %d;\n',max(nint));
        fprintf(fout,'\t-- Number of bits to represent memory address\n');
        fprintf(fout,'\tconstant N_BIT_ADDR_PWAS : integer := %d;\n',sum(nint));
        fprintf(fout,'\t-- Number of bits for the multipliers\n');
        fprintf(fout,'\tconstant N_BIT_MUL_PWAS : integer := %d;\n',nbit_mul);
        fprintf(fout,'\t\n');
        fprintf(fout,'\t-- Definition of new types\n');
        fprintf(fout,'\t\n');
        if ndim == 1
            if strcmpi(inputRepresentation,'signed')
                fprintf(fout,'\ttype x_matrix_pwas is array(1 downto 0) of signed(N_BIT_PWAS-1 downto 0);\n');
            else
                fprintf(fout,'\ttype x_matrix_pwas is array(1 downto 0) of signed(N_BIT_PWAS downto 0);\n');
            end
        else
            if strcmpi(inputRepresentation,'signed')
                fprintf(fout,'\ttype x_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_PWAS-1 downto 0);\n');
            else
                fprintf(fout,'\ttype x_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_PWAS downto 0);\n');
            end
        end
        fprintf(fout,'\ttype int_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of unsigned(N_INT_MAX_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype dec_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of unsigned(N_BIT_COEFF_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype y_matrix_pwas is array(N_FUN_PWAS-1 downto 0) of signed(N_BIT_OUT_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype mul_in_matrix_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_MUL_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype mul_out_matrix_pwas is array(N_DIM_PWAS downto 0) of signed(2*N_BIT_MUL_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype mu_matrix_pwas is array(N_DIM_PWAS downto 0) of unsigned(N_BIT_COEFF_PWAS-1 downto 0);\n');
        if ndim == 1
            fprintf(fout,'\ttype coeff_matrix_pwas is array(1 downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n');
        else
            fprintf(fout,'\ttype coeff_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n');
        end
        if ndim == 1
            fprintf(fout,'\ttype coeff_max_matrix_pwas is array(1 downto 0) of signed(2*N_BIT_MUL_PWAS-1 downto 0);\n');
        else
            fprintf(fout,'\ttype coeff_max_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(2*N_BIT_MUL_PWAS-1 downto 0);\n');
        end
        fprintf(fout,'\ttype d_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of std_logic_vector(N_DIM_PWAS-1 downto 0);\n');
        fprintf(fout,'\ttype addr_matrix_pwas is array(N_DIM_PWAS downto 0) of unsigned(N_BIT_ADDR_PWAS-1 downto 0);\n');
        if strcmpi(outputRepresentation,'signed')
            fprintf(fout,'\ttype mem_element_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n');
        else
            fprintf(fout,'\ttype mem_element_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_COEFF_PWAS downto 0);\n');
        end
        fprintf(fout,'\ttype mem_matrix_pwas is array(N_FUN_PWAS-1 downto 0) of mem_element_pwas;\n');
        if ndim == 1
            fprintf(fout,'\ttype integer_matrix_pwas is array(1 downto 0) of integer;\n');
        else
            fprintf(fout,'\ttype integer_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of integer;\n');
        end
        fprintf(fout,'\ttype integer_matrix_out_pwas is array(N_FUN_PWAS-1 downto 0) of integer;\n');
        if ~uniform
            if strcmpi(inputRepresentation,'signed')
                for i = 1:ndim
                    fprintf(fout,'\ttype p%d_matrix_pwas is array(%d downto 0) of signed(N_BIT_PWAS-1 downto 0);\n',i-1,np(i)-1);
                end
            else
                for i = 1:ndim
                    fprintf(fout,'\ttype p%d_matrix_pwas is array(%d downto 0) of signed(N_BIT_PWAS downto 0);\n',i-1,np(i)-1);
                end
            end
            for i = 1:ndim
            	fprintf(fout,'\ttype scale_matrix%d_a_pwas is array(%d downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n',i-1,np(i)-1);
                fprintf(fout,'\ttype scale_matrix%d_b_pwas is array(%d downto 0) of signed(N_BIT_PWAS-1 downto 0);\n',i-1,np(i)-1);
            end
        end
        
        fprintf(fout,'\n');
        fprintf(fout,'\t-- Declare other constants\n');
        fprintf(fout,'\n');
        fprintf(fout,'\t-- Number of bits to represent the integer part of the inputs after scaling\n');
        if ndim == 1
            fprintf(fout,'\tconstant N_BIT_INT_PWAS: integer_matrix_pwas := (0,');
            fprintf(fout,'%d);\n',nint(1));
        else
            fprintf(fout,'\tconstant N_BIT_INT_PWAS: integer_matrix_pwas := (');
            for i = ndim:-1:2
                fprintf(fout,'%d,',nint(i));
            end
            fprintf(fout,'%d);\n',nint(1));
        end
        if ndim == 1
            fprintf(fout,'\tconstant NP_PWAS: coeff_max_matrix_pwas := ((others => ''0''),');
            fprintf(fout,'"%s");\n',npbin(1).bin); 
        else
            fprintf(fout,'\tconstant NP_PWAS: coeff_max_matrix_pwas := (');
            for i = ndim:-1:2
                fprintf(fout,'"%s",',npbin(i).bin);
            end
            fprintf(fout,'"%s");\n',npbin(1).bin);
        end
        fprintf(fout,'\tconstant ZERO_PWAS : signed(2*N_BIT_MUL_PWAS-1 downto 0) := (others => ''0'');\n');
        fprintf(fout,'\t-- Point position of scaled terms A(x-b)\n');
        if ndim == 1
            fprintf(fout,'\tconstant POINT_POSITION_PWAS : integer_matrix_pwas := (0,');
            fprintf(fout,'%d);\n',App(1));
        else
            fprintf(fout,'\tconstant POINT_POSITION_PWAS : integer_matrix_pwas := (');
            for i = ndim:-1:2
                fprintf(fout,'%d,',App(i));
            end
            fprintf(fout,'%d);\n',App(1));
        end
        if uniform
            if ndim == 1
                fprintf(fout,'\tconstant scaleInput_A_PWAS : coeff_matrix_pwas := ((others => ''0''),');
                fprintf(fout,'"%s");\n',A(1).bin);     
            else
                fprintf(fout,'\tconstant scaleInput_A_PWAS : coeff_matrix_pwas := (');
                for i = ndim:-1:2
                    fprintf(fout,'"%s",',A(i).bin);
                end
                fprintf(fout,'"%s");\n',A(1).bin);
            end
            if ndim == 1
                fprintf(fout,'\tconstant scaleInput_b_PWAS : x_matrix_pwas := ((others => ''0''),');
                fprintf(fout,'"%s");\n',b(1).bin);
            else
                fprintf(fout,'\tconstant scaleInput_b_PWAS : x_matrix_pwas := (');
                for i = ndim:-1:2
                    fprintf(fout,'"%s",',b(i).bin);
                end
                fprintf(fout,'"%s");\n',b(1).bin);
            end
        else
            if ndim == 1
                fprintf(fout,'\tconstant NUMPART_PWAS: integer_matrix_pwas := (0,');
                fprintf(fout,'%d);\n',np(1));        
            else
                fprintf(fout,'\tconstant NUMPART_PWAS: integer_matrix_pwas := (');
                for i = ndim:-1:2
                    fprintf(fout,'%d,',np(i));
                end
                fprintf(fout,'%d);\n',np(1));
            end
            for i = 1:ndim
                fprintf(fout,'\tconstant scaleInput%d_A_PWAS : scale_matrix%d_a_pwas := (',i-1,i-1);
                for j = np(i):-1:2
                    fprintf(fout,'"%s",',A{i}(j).bin);
                end
                fprintf(fout,'"%s");\n',A{i}(1).bin);
            end
            for i = 1:ndim
                fprintf(fout,'\tconstant scaleInput%d_b_PWAS : scale_matrix%d_b_pwas := (',i-1,i-1);
                for j = np(i):-1:2
                    fprintf(fout,'"%s",',b{i}(j).bin);
                end
                fprintf(fout,'"%s");\n',b{i}(1).bin);
            end
            
        end
        fprintf(fout,'\n');
        if nfun > 1
            fprintf(fout,'\tconstant START_BIT_PWAS : integer_matrix_out_pwas := (');
        else
            fprintf(fout,'\tconstant START_BIT_PWAS : integer := (');
        end
        for i = nfun:-1:2
            fprintf(fout,'%d,',fpp(i)+nbit_coeff-1);
        end
        fprintf(fout,'%d);\n',fpp(1)+nbit_coeff-1);
        fprintf(fout,'\n');
        
        if ~uniform
            for i = 1:ndim
                fprintf(fout,'\tconstant P%d_PWAS : p%d_matrix_pwas := (',i-1,i-1);
                for j = np(i)+1:-1:3
                    fprintf(fout,'"%s",',Pbin{i}(j).bin);
                end
                fprintf(fout,'"%s");\n',Pbin{i}(2).bin);
            end
        end
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'pwasFunctionPackage.vhd']);