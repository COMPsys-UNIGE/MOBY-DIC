%% WRITE MEMORY
% Writes the VHDL code implementing block Memory in pwag circuit
%
% SYNTAX
%
% writeMemory(architecture,memorybin,folder)
%
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

function writeMemory(object,Hint,Kint,Fint,Gint,nbitMem,nbitcoeff,folder)

M = [Hint Kint; Fint Gint];
rowSize = size(M,2);
columnSize = size(M,1);


fin = fopen([getvhdlpath(),'pwagFunction/Memory.vhd'],'r');
fout = fopen([folder,'Memory.vhd'],'w');

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
        %             for i = 1:size(memorybin{1},1)
        %                 fprintf(fout,'\t\t\t\twhen "');
        %                 fprintf(fout,'%s',dec2bin(i-1,bit_addr));
        %                 fprintf(fout,'" => res <= "%s";\n',memorybin{1}(i,:));
        %             end
        
        %%%%% NEW BY TEO
        % Memory is organized as [H K;
        %                         F G];
        
        fprintf(fout,'type ROM_row is array (0 to %d) of signed(N_BIT_COEFF_PWAG-1 downto 0);\n\n',rowSize-1);
        fprintf(fout,'type ROM_type is array (0 to %d) of ROM_row;\n\n',columnSize-1);
        fprintf(fout,'constant ROM : ROM_type:=(\n');
        
        for i = 1:columnSize-1
            fprintf(fout,'(');
            for j=1:rowSize-1
                numi = decimal2signed(M(i,j),nbitcoeff,0);
                fprintf(fout,'"%s", ',numi.bin);
            end
            numi = decimal2signed(M(i,end),nbitcoeff,0);
            fprintf(fout,'"%s"),\n ',numi.bin);
        end
        fprintf(fout,'(');
        for j=1:rowSize-1
            numi = decimal2signed(M(end,j),nbitcoeff,0);
            fprintf(fout,'"%s", ',numi.bin);
        end
        numi = decimal2signed(M(end,end),nbitcoeff,0);
        fprintf(fout,'"%s"));\n ',numi.bin);
        
        
        %%%%%
        
        
        
        fprintf(fout,'\n');
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'Memory.vhd']);

end