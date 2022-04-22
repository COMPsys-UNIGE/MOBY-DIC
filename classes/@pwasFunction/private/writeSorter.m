function writeSorter(object,options)
% writeSorter
% Writes the VHDL code implementing block Sorter 
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
folder = options.folder;

if ndim > 4
    fin = fopen([getvhdlpath(),'pwasFunction/Sorter8.vhd'],'r');
    ndimmax = 8;
elseif ndim > 2
    fin = fopen([getvhdlpath(),'pwasFunction/Sorter4.vhd'],'r');
    ndimmax = 4;
else
    fin = fopen([getvhdlpath(),'pwasFunction/Sorter2.vhd'],'r');
    ndimmax = 2;
end

fout = fopen([folder,'Sorter.vhd'],'w');

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
        
        for i = 1:ndim
            fprintf(fout,'\tsignal_x%d <= x(%d);\n',i,i-1);
        end
        for i = ndim+1:ndimmax
            fprintf(fout,'\tsignal_x%d <= (others => ''0'');\n',i);
        end
        
        for i = 1:ndim
            fprintf(fout,'\ty(%d) <= signal_y%d;\n',i-1,i);
        end
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'Scale.vhd']);
