function writeMemory(object,idx,bits,opts,scaleOutput)
% writeMemory
% Writes the VHDL code implementing block Memory in pwasFunction circuit
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

% Extract input arguments
folder = opts.folder;
nbit_coeff = opts.coeffResolution;
outputRange = opts.outputRange;
outputRepresentation = opts.outputRepresentation;
nfun = numel(idx);
ndim = object.getDomainDimensions();
np = object.getNumberOfPartitions();
bits_int = bits.nint;

% Weights
w = object.getWeights();
w = w(:,idx);

wbin = cell(1,nfun);

% Loop oncodomain dimensions
for id = 1:numel(idx)
    
    % Normalize weights 
    wcur = scaleOutput.A(id)*w(:,id)+scaleOutput.b(id);
    
    % If the number of partitions is not the maximum allowed by the number
    % of bits of the integer part, zeros are inserted in
    % correspondance of the missing weights
    if prod(np+1) < prod(2.^bits_int)
        
        % Useful to avoid that the memory is filled with zeros in cells
        % which will be never accessed
        dims = 2.^bits_int;
        dims(end) = np(end)+1;
        
        if numel(np) == 1
            fillMatrix = zeros(dims,1);
        else
            fillMatrix = zeros(dims');
        end
        
        ind = cell(ndim,1);
        
        for n = 1:ndim
            ind{n} = (1:np(n)+1);
        end
        
        if numel(np) == 1
            fillMatrix(ind{:}) = full(wcur)';
        else
            fillMatrix(ind{:}) = reshape(full(wcur),np'+1);
        end
        wcur = fillMatrix(:);
        
    end
    
    % Output range
    fmin = outputRange.min(id);
    fmax = outputRange.max(id);
    if strcmpi(outputRepresentation,'signed')
        fminb = decimal2signed(fmin,nbit_coeff);
        fmaxb = decimal2signed(fmax,nbit_coeff);
    else
        fminb = decimal2unsigned(fmin,nbit_coeff);
        fmaxb = decimal2unsigned(fmax,nbit_coeff);
    end
    fpp = min(fminb.pointposition,fmaxb.pointposition);
    
    if strcmpi(outputRepresentation,'signed')
        wbin{id} = decimal2signed(wcur,nbit_coeff,fpp);
    else
        wbin{id} = decimal2unsigned(wcur,nbit_coeff,fpp);
    end
end

if strcmpi(opts.architecture,'fast')
    
    
    fin = fopen([getvhdlpath(),'pwasFunction/Memory.vhd'],'r');
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
            
            
            
            fprintf(fout,'\ttype ROM_type is array (0 to %d) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n\n',numel(wbin{1})-1);
            
            for id = 1:nfun
                fprintf(fout,'\tconstant ROM_%d : ROM_type:=(\n\t',id-1);
                for i = 1:numel(wbin{id})-1
                    fprintf(fout,'"%s", ',wbin{id}(i).bin);
                    if mod(i,10) == 0
                        fprintf(fout,'\n\t');
                    end
                end
                fprintf(fout,'"%s");\n',wbin{id}(end).bin);
                fprintf(fout,'\n');
            end
            
            fprintf(fout,'begin\n');
            fprintf(fout,'\n');
            
            fprintf(fout,'\tproc_mem: process(clk,reset)\n');
            fprintf(fout,'\tbegin\n');
            fprintf(fout,'\t\tif reset = ''0'' then\n');
            fprintf(fout,'\t\t\tfor i in 0 to N_FUN_PWAS-1 loop\n');
            fprintf(fout,'\t\t\t\tfor j in 0 to N_DIM_PWAS loop\n');
            fprintf(fout,'\t\t\t\t\tres(i)(j) <= (others => ''0'');\n');
            fprintf(fout,'\t\t\t\tend loop;\n');
            fprintf(fout,'\t\t\tend loop;\n');
            fprintf(fout,'\t\telsif rising_edge(clk) and start = ''1'' then\n');
            fprintf(fout,'\t\t\tfor j in 0 to N_DIM_PWAS loop\n');
            if strcmpi(outputRepresentation,'signed')
                for i = 1:nfun
                    fprintf(fout,'\t\t\t\tres(%d)(j) <= ROM_%d(to_integer(addr(j)));\n',i-1,i-1);
                end
            else
                for i = 1:nfun
                    fprintf(fout,'\t\t\t\tres(%d)(j) <= ''0''&ROM_%d(to_integer(addr(j)));\n',i-1,i-1);
                end
            end
            fprintf(fout,'\t\t\tend loop;\n');
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
    
    disp(['Generated file ',folder,'Memory.vhd']);
    
else
    
    fin = fopen([getvhdlpath(),'pwasFunction/Memory_ser.vhd'],'r');
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
            
            fprintf(fout,'\tsignal state : integer range 0 to N_DIM_PWAS+1;\n');
            
            fprintf(fout,'\ttype ROM_type is array (0 to %d) of signed(N_BIT_COEFF_PWAS-1 downto 0);\n\n',numel(wbin{1})-1);
            
            for id = 1:nfun
                fprintf(fout,'\tconstant ROM_%d : ROM_type:=(\n\t',id-1);
                for i = 1:numel(wbin{id})-1
                    fprintf(fout,'"%s", ',wbin{id}(i).bin);
                    if mod(i,10) == 0
                        fprintf(fout,'\n\t');
                    end
                end
                fprintf(fout,'"%s");\n',wbin{id}(end).bin);
                fprintf(fout,'\n');
            end
            
            fprintf(fout,'begin\n');
            fprintf(fout,'\n');
            fprintf(fout,'\tproc_mem: process(clk,reset,state,start)\n');
            fprintf(fout,'\tbegin\n');
            fprintf(fout,'\t\tif reset = ''0'' then\n');
            fprintf(fout,'\t\t\tfor i in 0 to N_FUN_PWAS-1 loop\n');
            fprintf(fout,'\t\t\t\tfor j in 0 to N_DIM_PWAS loop\n');
            fprintf(fout,'\t\t\t\t\tres(i)(j) <= (others => ''0'');\n');
            fprintf(fout,'\t\t\t\tend loop;\n');
            fprintf(fout,'\t\t\tend loop;	\n');
            fprintf(fout,'\t\telsif rising_edge(clk) then\n');
            fprintf(fout,'\t\t\tif (state = 0 and start = ''1'') or (state > 0 and state < N_DIM_PWAS+1) then\n');
            if strcmpi(outputRepresentation,'signed')
                for i = 1:nfun
                    fprintf(fout,'\t\t\t\tres(%d)(state) <= ROM_%d(to_integer(addr(state)));\n',i-1,i-1);
                end
            else
                for i = 1:nfun
                    fprintf(fout,'\t\t\t\tres(%d)(state) <= ''0''&ROM_%d(to_integer(addr(state)));\n',i-1,i-1);
                end
            end
            fprintf(fout,'\t\t\tend if;\n');
            fprintf(fout,'\t\tend if;\n');
            fprintf(fout,'\tend process;\n');
            
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

end