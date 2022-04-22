function writeScale(object,optout,bits)
% writeScale
% Writes the VHDL code implementing block Scale
%
% SYNTAX
% This is a private method

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Matteo Lodi (matteo.lodi@edu.unige.it)
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
% Boston, MA  02111-1307  USA

ndim = object.getDomainDimensions();
folder = optout.folder;
np = object.getNumberOfPartitions();

nint = bits.nint;
nintmax = max(nint);

fin = fopen([getvhdlpath(),'pwasFunction/Scale_nu.vhd'],'r');
fout = fopen([folder,'Scale.vhd'],'w');

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
        
        fprintf(fout,'\t-- In this process the indexes of the hyper-rectangle containing the point are found\n');
        fprintf(fout,'\tproc_find_index: process(reset,x)\n');
        fprintf(fout,'\tbegin\n');
        
        fprintf(fout,'\t\tif reset = ''0'' then\n');
        fprintf(fout,'\t\t\tfor i in 0 to N_DIM_PWAS-1 loop\n');
        fprintf(fout,'\t\t\t\tsignal_scaleInput_A_PWAS(i) <= (others => ''0'');\n');
        fprintf(fout,'\t\t\t\tsignal_scaleInput_b_PWAS(i) <= (others => ''0'');\n');
        fprintf(fout,'\t\t\t\tint(i) <= (others => ''0'');\n');
        fprintf(fout,'\t\t\tend loop;\n');
        fprintf(fout,'\t\telse\n');
		
        for i = 1:ndim
            fprintf(fout,'\t\t\tif x(%d) < P0_PWAS(0) then\n',i-1);
            fprintf(fout,'\t\t\t\tsignal_scaleInput_A_PWAS(%d) <= scaleInput%d_A_PWAS(0);\n',i-1,i-1);
			fprintf(fout,'\t\t\t\tsignal_scaleInput_b_PWAS(%d) <= scaleInput%d_b_PWAS(0);\n',i-1,i-1);
			int = decimal2unsigned(0,nintmax,0);
            fprintf(fout,'\t\t\t\tint(%d) <= "%s";\n',i-1,int.bin);
            
            for j = 1:np(i)-2
                fprintf(fout,'\t\t\telsif x(%d) < P0_PWAS(%d) then\n',i-1,j);
                fprintf(fout,'\t\t\t\tsignal_scaleInput_A_PWAS(%d) <= scaleInput%d_A_PWAS(%d);\n',i-1,i-1,j);
                fprintf(fout,'\t\t\t\tsignal_scaleInput_b_PWAS(%d) <= scaleInput%d_b_PWAS(%d);\n',i-1,i-1,j);
                int = decimal2unsigned(j,nintmax,0);
                fprintf(fout,'\t\t\t\tint(%d) <= "%s";\n',i-1,int.bin);
            end
            
            fprintf(fout,'\t\t\telse\n');
			fprintf(fout,'\t\t\t\tsignal_scaleInput_A_PWAS(%d) <= scaleInput%d_A_PWAS(%d);\n',i-1,i-1,np(i)-1);
			fprintf(fout,'\t\t\t\tsignal_scaleInput_b_PWAS(%d) <= scaleInput%d_b_PWAS(%d);\n',i-1,i-1,np(i)-1);
            int = decimal2unsigned(np(i)-1,nintmax,0);
			fprintf(fout,'\t\t\t\tint(%d) <= "%s";\n',i-1,int.bin);
            fprintf(fout,'\t\t\tend if;\n');
        end

	fprintf(fout,'\t\tend if;\n');
	
	fprintf(fout,'\tend process;\n');
            
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'Scale.vhd']);
