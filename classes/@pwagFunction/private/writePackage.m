function writePackage(object,nbit,nbit_coeff,nbitout,nbitAddr,nbitState,nintcoeff,alpha,beta,idx,folder)
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
% folder = opts.folder;
% nbit = bits.nbit;
% nbit_coeff = bits.nbit_coeff;
% nbitout = bits.nbitout;

% Number of domain and codomain dimensions
ndim = object.getDomainDimensions();
nfun = object.getCodomainDimensions();


% Actual number of bits to represent input: if input representation is
% unsigned, one more bit is necessary because the multipliers are signed
% and therefore a '0' must be added as MSB of the input
nbit_eff = nbit;

% Idem for the output
nbit_coeff_eff = nbit_coeff;


% Number of bits for the multiplier
nbit_mul = max(nbit_eff,nbit_coeff_eff);


% Open files
fin = fopen([getvhdlpath(),'pwagFunction/pwagFunctionPackage.vhd'],'r');
fout = fopen([folder,'pwagFunctionPackage.vhd'],'w');

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
        fprintf(fout,'\tconstant N_DIM_PWAG : integer := %d;\n',ndim);
        fprintf(fout,'\t-- Number of codomain dimensions\n');
        fprintf(fout,'\tconstant N_FUN_PWAG : integer := %d;\n',numel(idx));       
        fprintf(fout,'\t-- Number of bits to represent the inputs\n');
        fprintf(fout,'\tconstant N_BIT_PWAG : integer := %d;\n',nbit);
        fprintf(fout,'\t-- Number of bits to represent the outputs\n');
        fprintf(fout,'\tconstant N_BIT_OUT_PWAG : integer := %d;\n',nbitout);
        fprintf(fout,'\t-- Number of bits to represent the coefficients\n');
        fprintf(fout,'\tconstant N_BIT_COEFF_PWAG : integer := %d;\n',nbit_coeff);
        fprintf(fout,'\t-- Number of bits to represent memory address\n');
        fprintf(fout,'\tconstant N_BIT_ADDR_PWAG : integer := %d;\n',nbitAddr);
        fprintf(fout,'\t-- Number of bits for the multipliers\n');
        fprintf(fout,'\tconstant N_BIT_MUL_PWAG : integer := %d;\n',nbit_mul);
        fprintf(fout,'\t-- Number of bits for the stateReg\n');
        fprintf(fout,'\tconstant N_BIT_STATE_PWAG : integer := %d;\n',nbitState);
        fprintf(fout,'\t\n');
        fprintf(fout,'\t-- Definition of new types\n');
        fprintf(fout,'\t\n');
%         if strcmpi(inputRepresentation,'signed')
            fprintf(fout,'\ttype x_matrix_PWAG is array(N_DIM_PWAG-1 downto 0) of signed(N_BIT_PWAG-1 downto 0);\n');
%         else
%             fprintf(fout,'\ttype x_matrix is array(N_DIM-1 downto 0) of unsigned(N_BIT downto 0);\n');
%         end
        fprintf(fout,'\ttype y_matrix_PWAG is array(N_FUN_PWAG-1 downto 0) of signed(N_BIT_OUT_PWAG-1 downto 0);\n');
        fprintf(fout,'\ttype mul_in_matrix_PWAG is array(N_DIM_PWAG-1 downto 0) of signed(N_BIT_MUL_PWAG-1 downto 0);\n');
        fprintf(fout,'\ttype mul_out_matrix_PWAG is array(N_DIM_PWAG-1 downto 0) of signed(2*N_BIT_MUL_PWAG-1 downto 0);\n');
        fprintf(fout,'\ttype addr_array_PWAG is array(N_FUN_PWAG-1 downto 0) of unsigned(N_BIT_ADDR_PWAG-1 downto 0);\n');
        fprintf(fout,'\ttype memory_out_matrix_PWAG is array(N_DIM_PWAG downto 0) of signed(N_BIT_COEFF_PWAG-1 downto 0);\n');
        fprintf(fout,'\ttype int_vect_PWAG is array(0 to N_DIM_PWAG) of integer;\n');
        fprintf(fout,'\tconstant preshift : int_vect_PWAG := (');
        nshiftMin = min(nintcoeff);
        maxShift = max(nintcoeff-nshiftMin);
        for i=1:numel(nintcoeff)-1
            fprintf(fout,'%d,',nintcoeff(i)-nshiftMin);
        end
        ndeccoeff = nbit_coeff-nintcoeff;
        fprintf(fout,'%d);\n',nintcoeff(end)-nshiftMin);
        fprintf(fout,'\tconstant N_BIT_ADDER_PWAG : integer := %d;\n',nbit+nbit_coeff+ceil(log2(ndim))+maxShift-nshiftMin);
        fprintf(fout,'\tconstant N_DEC_SUM_OUT_PWAG : integer := %d;\n',max(ndeccoeff));
        alphabin = decimal2signed(alpha,nbit_coeff);
        fprintf(fout,'\t constant alpha_PWAG : signed(N_BIT_COEFF_PWAG-1 downto 0) := "%s";\n',alphabin.bin);
         betabin = decimal2signed(beta,nbit_coeff+nbit+nbit_coeff+ceil(log2(ndim))+maxShift-nshiftMin,max(ndeccoeff) + alphabin.pointposition);
        fprintf(fout,'\t constant beta_PWAG : signed(N_BIT_COEFF_PWAG+N_BIT_ADDER_PWAG-1 downto 0) := "%s";\n',betabin.bin);
        satPos = decimal2signed(2^(nbitout-1)-1,nbit_coeff+nbit+nbit_coeff+ceil(log2(ndim))+maxShift-nshiftMin,max(ndeccoeff) + alphabin.pointposition); 
        fprintf(fout,'\t constant maxVal_PWAG : signed(N_BIT_COEFF_PWAG+N_BIT_ADDER_PWAG-1 downto 0) := "%s";\n',satPos.bin);
        satNeg = decimal2signed(-(2^(nbitout-1)),nbit_coeff+nbit+nbit_coeff+ceil(log2(ndim))+maxShift-nshiftMin,max(ndeccoeff) + alphabin.pointposition); 
        fprintf(fout,'\t constant minVal_PWAG : signed(N_BIT_COEFF_PWAG+N_BIT_ADDER_PWAG-1 downto 0) := "%s";\n',satNeg.bin);

         fprintf(fout,'\t constant rescaleOutDecPosistion_PWAG : integer := %d;',alphabin.pointposition + max(ndeccoeff));
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'pwagFunctionPackage.vhd']);