function writePackage(object,opts,perf)
% writePackage
% Writes the VHDL code implementing package vsPackage
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

nbit = opts.inputResolution;
nbitout = opts.outputResolution;
ndim = object.getNumberOfDimensions();
nbitmul = perf.nbitmul;
folder = opts.folder;
frequency = opts.frequency;
Ts = object.getSamplingTime();

mu = object.getInputTimeWindow();
my = object.getOutputTimeWindow();
mz = object.getAutoregressiveTimeWindow();

defz = opts.initialCondition;

if strcmpi(opts.inputRepresentation,'signed')
    defz = decimal2signed(defz,nbit,0);
else
    defz = decimal2unsigned(defz,nbit,0);
end

clkcycles = Ts*frequency;

fun = object.getFunction();
npwas = numel(fun);

% Open files
fin = fopen([getvhdlpath(),'virtualSensor/vsPackage.vhd'],'r');
fout = fopen([folder,'vsPackage.vhd'],'w');

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
        
        fprintf(fout,'\t-- Number of bits to represent the inputs\n');
        fprintf(fout,'\tconstant N_BIT_VS : integer := %d;\n',nbit);
        fprintf(fout,'\t-- Number of bits to represent the outputs\n');
        fprintf(fout,'\tconstant N_BIT_OUT_VS : integer := %d;\n',nbitout);
        
        fprintf(fout,'\tconstant N_DIM_VS : integer := %d;\n',ndim);
        fprintf(fout,'\tconstant N_BIT_MUL_VS : integer := %d;\n',nbitmul);
        fprintf(fout,'\tconstant N_FUN_VS : integer := 1;\n');
        if mz > 0
            fprintf(fout,'\tconstant DEFAULT_Z : std_logic_vector(N_BIT_VS-1 downto 0) := "%s";\n',defz.bin);
        end
        fprintf(fout,'\tconstant sampling_latency_VS : integer := %d;\n',clkcycles);
        
        if object.isCurrent
            if mu > 0
                fprintf(fout,'\tconstant MU_VS : integer := %d;\n',max(mu));
            else
                fprintf(fout,'\tconstant MU_VS : integer := 0;\n');
            end
            if my > 0
                fprintf(fout,'\tconstant MY_VS : integer := %d;\n',max(my));
            end
            fprintf(fout,'\tconstant MZ_VS : integer := %d;\n',mz);
            fprintf(fout,'\tconstant MMAX_VS : integer := %d;\n',max([max(mu(:));max(my(:))]));
        else
            if mu > 0
                fprintf(fout,'\tconstant MU_VS : integer := %d;\n',max(mu)+1);
            else
                fprintf(fout,'\tconstant MU_VS : integer := 0;\n');
            end
            if my > 0
                fprintf(fout,'\tconstant MY_VS : integer := %d;\n',max(my)+1);
            end
            fprintf(fout,'\tconstant MZ_VS : integer := %d;\n',mz);
            fprintf(fout,'\tconstant MMAX_VS : integer := %d;\n',max([max(mu(:)+1);max(my(:)+1)]));
        end
        
        fprintf(fout,'\t-- Definition of new types\n');
        
        fprintf(fout,'\ttype buf_u_matrix_vs is array(MU_VS-1 downto 0) of std_logic_vector(N_BIT_VS-1 downto 0);\n');
        fprintf(fout,'\ttype buf_y_matrix_vs is array(MY_VS-1 downto 0) of std_logic_vector(N_BIT_VS-1 downto 0);\n');
        fprintf(fout,'\ttype buf_z_matrix_vs is array(MZ_VS-1 downto 0) of std_logic_vector(N_BIT_VS-1 downto 0);\n');
        
        for i = 1:npwas
            ndim = fun(i).getDomainDimensions();
            if ndim == 1
                if strcmpi(opts.inputRepresentation,'signed')
                    fprintf(fout,'\ttype x_matrix%d_vs is array(1 downto 0) of signed(N_BIT_VS-1 downto 0);\n',i);
                else
                    fprintf(fout,'\ttype x_matrix%d_vs is array(1 downto 0) of signed(N_BIT_VS downto 0);\n',i);
                end
            else
                if strcmpi(opts.inputRepresentation,'signed')
                    fprintf(fout,'\ttype x_matrix%d_vs is array(%d downto 0) of signed(N_BIT_VS-1 downto 0);\n',i,ndim-1);
                else
                    fprintf(fout,'\ttype x_matrix%d_vs is array(%d downto 0) of signed(N_BIT_VS downto 0);\n',i,ndim-1);
                end
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