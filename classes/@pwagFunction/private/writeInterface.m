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
inputRepresentation = opts.inputRepresentation;
ndim = object.getDomainDimensions();
nfun = object.getCodomainDimensions();

fin = fopen([getvhdlpath(),'pwagFunction/pwagFunction.vhd'],'r');
fout = fopen([folder,'pwagFunction.vhd'],'w');
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
            
            fprintf(fout,'entity pwagFunction is\n');
            fprintf(fout,'\tPort ( clk : in std_logic;\n');
            fprintf(fout,'\t        ce : in std_logic;\n');
            fprintf(fout,'\t     reset : in std_logic;\n');
            fprintf(fout,'\t     start : in std_logic;\n');
            for i = 1:ndim
                fprintf(fout,'\t       x%d : in std_logic_vector(N_BIT_PWAG-1 downto 0);\n',i);
            end
            for i = 1:nfun
                fprintf(fout,'\t       y%d : out std_logic_vector(N_BIT_OUT_PWAG-1 downto 0);\n',i);
            end
            fprintf(fout,'\t      done : out std_logic);\n');
            fprintf(fout,'end pwagFunction;\n');
            
        end
        
        if count_matlabgen == 2

            for i = 1:nfun
                fprintf(fout,'\ty%d <= std_logic_vector(signal_y_scale(%d));\n',i,i-1);
            end
            
                for i = 1:ndim
                    fprintf(fout,'\tsignal_x(%d) <= signed(x%d);\n',i-1,i);
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

disp(['Generated file ',folder,'pwagFunction.vhd']);
