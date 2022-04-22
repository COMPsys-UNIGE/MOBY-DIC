function writeInterface(object,opts)
% writeInterface
% Writes the VHDL code implementing block pwasFunction
%
% SYNTAX
% This is a private method
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Alberto Oliveri (alberto.oliveri@unige.it)
%
% Copyright is with:
%
% * Copyright (C) 2011-2011 University of Genoa, Italy.

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

folder = opts.folder;
nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nxref = numel(object.getTrackingVariable);
nu = object.getNumberOfInputs;

fin = fopen([getvhdlpath(),'ApproxMPCctrl/controllerInterface.vhd'],'r');
fout = fopen([folder,'controllerInterface.vhd'],'w');
count_matlabgen = 1;

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
        fprintf(fout,'\r\n');
        fprintf(fout,'\n');
        
        if count_matlabgen == 1
            
            fprintf(fout,'entity controllerInterface is\n');
            fprintf(fout,'\tPort ( clk : in std_logic;\n');
            fprintf(fout,'\t        ce : in std_logic;\n');
            fprintf(fout,'\t     reset : in std_logic;\n');
            
            for i = 1:nx
                fprintf(fout,'\t       x%d : in std_logic_vector(N_BIT_CTRL-1 downto 0);\n',i);
            end
            for i = 1:np
                fprintf(fout,'\t       p%d : in std_logic_vector(N_BIT_CTRL-1 downto 0);\n',i);
            end
            for i = 1:nd
                fprintf(fout,'\t       d%d : in std_logic_vector(N_BIT_CTRL-1 downto 0);\n',i);
            end
            for i = 1:nxref
                fprintf(fout,'\t       xref%d : in std_logic_vector(N_BIT_CTRL-1 downto 0);\n',i);
            end
            for i = 1:nu
                fprintf(fout,'\t       u%d : out std_logic_vector(N_BIT_OUT_CTRL-1 downto 0);\n',i);
            end
            fprintf(fout,'\t     sample : out std_logic;\n');
            fprintf(fout,'\t      done : out std_logic);\n');
            fprintf(fout,'end controllerInterface;\n');
            
        end
        
        if count_matlabgen == 2
            
            for i = 1:nu
                fprintf(fout,'\tu%d <= std_logic_vector(u_reg(%d));\n',i,i-1);
            end
            
            for i = 1:nx
                fprintf(fout,'\tsignal_x(%d) <= signed(x%d);\n',i-1,i);
            end
            for i = 1:np
                fprintf(fout,'\tsignal_x(%d+nx_CTRL) <= signed(p%d);\n',i-1,i);
            end
            for i = 1:nd
                fprintf(fout,'\tsignal_x(%d+nx_CTRL+npar_CTRL) <= signed(d%d);\n',i-1,i);
            end
            for i = 1:nxref
                fprintf(fout,'\tsignal_x(%d+nx_CTRL+npar_CTRL+nd_CTRL) <= signed(xref%d);\n',i-1,i);
            end
            
            
            
        end
        
        fprintf(fout,'\n');
        
        count_matlabgen = count_matlabgen + 1;
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'controllerInterface.vhd']);
