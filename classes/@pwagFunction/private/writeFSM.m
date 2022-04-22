%% WRITE FSM
% Writes the VHDL code implementing block FSM in pwag circuit
%
% SYNTAX
%
% writeFSM(architecture,state,memsize,nfun,folder)
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

function writeFSM(object,state,nbitmem,nfun,folder)

% type = architecture.type;

% if strcmpi(type,'serial')
    fin = fopen([getvhdlpath(),'pwagFunction/FSM.vhd'],'r');
    fout = fopen([folder,'FSM.vhd'],'w');
% else
%     fin = fopen([getvhdlpath(),'pwag/pwag_par/FSM.vhd'],'r');
%     fout = fopen([folder,'FSM.vhd'],'w');
% end


max_state = numel(state)-1;
bit_state = ceil(log2(max_state+1));

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
            
            % This way, a ROM is synthesized (slower)
            %                 for i=1:numel(state)
            %                     if ~state(i).leaf
            %                         fprintf(fout,'\t\t\t\twhen "%s" => Next_State <= "%s" + ("%s" & decision);\n',dec2bin(i-1,bit_state),...
            %                             dec2bin(state(i).next(1)-1,bit_state),repmat('0',1,bit_state-1));
            %                     end
            %                 end
            
            % This way, a FSM is synthesized (faster)
            for i=1:numel(state)
                if ~state(i).leaf
                    numi = decimal2unsigned(i-1,bit_state,0);
                    fprintf(fout,'\t\t\t\twhen "%s" => \n',numi.bin);
                    fprintf(fout,'\t\t\t\t\tif decision = ''0'' then \n');
                    numi = decimal2unsigned(state(i).next(1)-1,bit_state,0);
                    fprintf(fout,'\t\t\t\t\t\tNext_State <= "%s";\n',numi.bin);
                    fprintf(fout,'\t\t\t\t\telse\n');
                    numi = decimal2unsigned(state(i).next(2)-1,bit_state,0);
                    fprintf(fout,'\t\t\t\t\t\tNext_State <= "%s";\n',numi.bin);
                    fprintf(fout,'\t\t\t\t\tend if;\n');
                end
            end
            
            
        end
        
        if count_matlabgen == 2
            
            for i = 1:numel(state)
                 numi = decimal2unsigned(i-1,bit_state,0);
                fprintf(fout,'\t\t\twhen "%s" =>\n',numi.bin);
                if state(i).leaf==0
                    numi = decimal2unsigned(state(i).addrHK-1,nbitmem,0);
                    fprintf(fout,'\t\t\t\tsignal_addr(0) <= "%s";\n',numi.bin);
                    for j = 2:nfun
                        fprintf(fout,'\t\t\t\tsignal_addr(%d) <= (others => ''0'');\n',j-1);
                    end
                else
                    for j = 1:nfun
                        numi = decimal2unsigned(state(i).addrFG(j)-1,nbitmem,0);
                        fprintf(fout,'\t\t\t\tsignal_addr(%d) <= "%s";\n',j-1,numi.bin);
                    end
                end
                fprintf(fout,'\t\t\t\tsignal_leaf <= ''%d'';\n',state(i).leaf);
            end
            
            fprintf(fout,'\t\t\twhen others => \n');
            for j = 1:nfun
                fprintf(fout,'\t\t\t\tsignal_addr(%d) <= (others => ''0'');\n',j-1);
            end
            fprintf(fout,'\t\t\t\tsignal_leaf <= ''0'';\n');
            
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

disp(['Generated file ',folder,'FSM.vhd']);