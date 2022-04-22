function writeInterface(object,opts)
% writeInterface
% Writes the VHDL code implementing block virtualSensor
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

fin = fopen([getvhdlpath(),'virtualSensor/virtualSensor.vhd'],'r');
fout = fopen([folder,'virtualSensor.vhd'],'w');

nu = object.getNumberOfInputs();
ny = object.getNumberOfMeasurableOutputs();
mu = object.getInputTimeWindow();
my = object.getOutputTimeWindow();
mz = object.getAutoregressiveTimeWindow();
reducedComplexity = object.isReducedComplexity();
fun = object.getFunction();
npwas = numel(fun);

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
        
        fprintf(fout,'entity virtualSensor is\n');
        fprintf(fout,'\tPort ( clk : in std_logic;\n');
        fprintf(fout,'\t     reset : in std_logic;\n');
        for i = 1:nu
            fprintf(fout,'\t        u%d : in std_logic_vector(N_BIT_VS-1 downto 0);\n',i);
        end
        for i = 1:ny
            fprintf(fout,'\t        y%d : in std_logic_vector(N_BIT_VS-1 downto 0);\n',i);
        end
        fprintf(fout,'\t         z : out std_logic_vector(N_BIT_OUT_VS-1 downto 0);\n');
        fprintf(fout,'\t      done : out std_logic);\n');
        fprintf(fout,'end virtualSensor;\n');
        fprintf(fout,'\n');
        
        fprintf(fout,'architecture Behavioral of virtualSensor is\n');
        
        for i = 1:npwas
            fprintf(fout,'\tcomponent pwasFunction_%d is\n',i);
            fprintf(fout,'\t\tPort ( clk : in std_logic;\n');
            fprintf(fout,'\t\t     reset : in std_logic;\n');
            fprintf(fout,'\t\t     start : in std_logic;\n');
            for j = 1:fun(i).getDomainDimensions
                fprintf(fout,'\t\t        x%d : in std_logic_vector(N_BIT_VS-1 downto 0);\n',j);
            end
            fprintf(fout,'\t\t        y1 : out std_logic_vector(N_BIT_OUT_VS-1 downto 0);\n');
            fprintf(fout,'\t\t      done : out std_logic);\n');
            fprintf(fout,'\tend component;\n');
            fprintf(fout,'\n');
        end
        
        if nu ~= 0
            fprintf(fout,'\tcomponent bufferU is\n');
            fprintf(fout,'\t\tPort (    clk : in std_logic;\n');
            fprintf(fout,'\t\t        reset : in std_logic;\n');
            fprintf(fout,'\t\t       sample : in std_logic;\n');
            fprintf(fout,'\t\t          uin : in std_logic_vector(N_BIT_VS-1 downto 0);\n');
            fprintf(fout,'\t\t         uout : out buf_u_matrix_vs);\n');
            fprintf(fout,'\tend component;\n');
            fprintf(fout,'\n');
        end
        
        if ny ~= 0
            fprintf(fout,'\tcomponent bufferY is\n');
            fprintf(fout,'\t\tPort (    clk : in std_logic;\n');
            fprintf(fout,'\t\t        reset : in std_logic;\n');
            fprintf(fout,'\t\t       sample : in std_logic;\n');
            fprintf(fout,'\t\t          yin : in std_logic_vector(N_BIT_VS-1 downto 0);\n');
            fprintf(fout,'\t\t         yout : out buf_y_matrix_vs);\n');
            fprintf(fout,'\tend component;\n');
            fprintf(fout,'\n');
        end
        
        if mz ~= 0
            fprintf(fout,'\tcomponent bufferZ is\n');
            fprintf(fout,'\t\tGeneric ( default : std_logic_vector(N_BIT_VS-1 downto 0) := (others => ''0''));\n');
            fprintf(fout,'\t\tPort (    clk : in std_logic;\n');
            fprintf(fout,'\t\t        reset : in std_logic;\n');
            fprintf(fout,'\t\t       sample : in std_logic;\n');
            fprintf(fout,'\t\t          zin : in std_logic_vector(N_BIT_VS-1 downto 0);\n');
            fprintf(fout,'\t\t         zout : out buf_z_matrix_vs);\n');
            fprintf(fout,'\tend component;\n');
            fprintf(fout,'\n');
        end
        
        for i = 1:npwas
            fprintf(fout,'\tsignal signal_done%d: std_logic;\n',i);
            fprintf(fout,'\tsignal signal_finish%d: std_logic;\n',i);
            fprintf(fout,'\tsignal signal_z%d : std_logic_vector(N_BIT_OUT_VS-1 downto 0);\n',i);
        end
        fprintf(fout,'\tsignal signal_finish: std_logic;\n');
        for i = 1:nu
            fprintf(fout,'\tsignal u%dpast : buf_u_matrix_vs;\n',i);
        end
        for i = 1:ny
            fprintf(fout,'\tsignal y%dpast : buf_y_matrix_vs;\n',i);
        end
        if mz ~= 0
            fprintf(fout,'\tsignal zpast : buf_z_matrix_vs;\n');
        end
        
        fprintf(fout,'\tsignal signal_start : std_logic;\n');
        fprintf(fout,'\tsignal signal_done : std_logic;\n');
        fprintf(fout,'\tsignal signal_sample : std_logic;\n');
        fprintf(fout,'\tsignal signal_full : std_logic;\n');
        fprintf(fout,'\tsignal signal_z : std_logic_vector(N_BIT_OUT_VS-1 downto 0);\n');
        fprintf(fout,'\tsignal counter : integer range 0 to sampling_latency_VS-1;\n');
        fprintf(fout,'\tsignal counter_full : integer range 0 to MMAX_VS;\n');
        fprintf(fout,'\n');
        fprintf(fout,'begin\n');
        
        for ii = 1:npwas
            fprintf(fout,'\tinst_pwasFunction_%d: pwasFunction_%d\n',ii,ii);
            fprintf(fout,'\t\tport map ( clk => clk,\n');
            fprintf(fout,'\t\t         reset => reset,\n');
            fprintf(fout,'\t\t         start => signal_start,\n');
            
            if ~reducedComplexity
                
                k = 1;
                for i = 1:nu
                    for j = 1:mu
                        if object.isCurrent
                            fprintf(fout,'\t\t            x%d => u%dpast(%d),\n',k,i,j-1);
                        else
                            fprintf(fout,'\t\t            x%d => u%dpast(%d),\n',k,i,j);
                        end
                        k = k+1;
                    end
                end
                for i = 1:ny
                    for j = 1:my
                        if object.isCurrent
                            fprintf(fout,'\t\t            x%d => y%dpast(%d),\n',k,i,j-1);
                        else
                            fprintf(fout,'\t\t            x%d => y%dpast(%d),\n',k,i,j);
                        end
                        k = k+1;
                    end
                end
                for j = 1:mz
                    fprintf(fout,'\t\t            x%d => zpast(%d),\n',k,j-1);
                    k = k+1;
                end
                
            else
                
                unames = cell(npwas,nu);
                ynames = cell(npwas,ny);
                znames = cell(npwas,1);
                
                if object.isCurrent
                    for i = 1:nu
                        for j = 1:mu(i)
                            unames{j,i} = ['u',num2str(i),'past(',num2str(j-1),')'];
                        end
                    end
                    for i = 1:ny
                        for j = 1:my(i)
                            ynames{j,i} = ['y',num2str(i),'past(',num2str(j-1),')'];
                        end
                    end
                    for j = 2:mz+1
                        znames{j,1} = ['zpast(',num2str(j-2),')'];
                    end
                else
                    for i = 1:nu
                        for j = 1:mu(i)
                            unames{j,i} = ['u',num2str(i),'past(',num2str(j),')'];
                        end
                    end
                    for i = 1:ny
                        for j = 1:my(i)
                            ynames{j,i} = ['y',num2str(i),'past(',num2str(j),')'];
                        end
                    end
                    for j = 1:mz
                        znames{j,1} = ['zpast(',num2str(j-1),')'];
                    end
                    
                end
                k = 1;
                for i = 1:nu
                    if ~isempty(unames{ii,1})
                        fprintf(fout,'\t\t            x%d => %s,\n',k,unames{ii,i});
                        k = k+1;
                    end
                end
                for i = 1:ny
                    if ~isempty(ynames{ii,i})
                        fprintf(fout,'\t\t            x%d => %s,\n',k,ynames{ii,i});
                        k = k+1;
                    end
                end
                for j = 1:mz
                    if ~isempty(znames{ii,1})
                        fprintf(fout,'\t\t            x%d => %s,\n',k,znames{ii,1});
                        k = k+1;
                    end
                end
            end
            
            fprintf(fout,'\t\t            y1 => signal_z%d,\n',ii);
            fprintf(fout,'\t\t          done => signal_done%d);\n',ii);
            fprintf(fout,'\n');
        end
        
        for i = 1:nu
            fprintf(fout,'\tinst_bufferU%d: bufferU\n',i);
            fprintf(fout,'\t\tport map ( clk => clk,\n');
            fprintf(fout,'\t\t         reset => reset,\n');
            fprintf(fout,'\t\t        sample => signal_sample,\n');
            fprintf(fout,'\t\t           uin => u%d,\n',i);
            fprintf(fout,'\t\t          uout => u%dpast);\n',i);
            fprintf(fout,'\n');
        end
        
        for i = 1:ny
            fprintf(fout,'\tinst_bufferY%d: bufferY\n',i);
            fprintf(fout,'\t\tport map ( clk => clk,\n');
            fprintf(fout,'\t\t         reset => reset,\n');
            fprintf(fout,'\t\t        sample => signal_sample,\n');
            fprintf(fout,'\t\t           yin => y%d,\n',i);
            fprintf(fout,'\t\t          yout => y%dpast);\n',i);
            fprintf(fout,'\n');
        end
        
        if mz > 0
            fprintf(fout,'\tinst_bufferZ: bufferZ\n');
            fprintf(fout,'\t\tgeneric map ( default => DEFAULT_Z)\n');
            fprintf(fout,'\t\tport map ( clk => clk,\n');
            fprintf(fout,'\t\t         reset => reset,\n');
            fprintf(fout,'\t\t        sample => signal_done,\n');
            fprintf(fout,'\t\t           zin => signal_z,\n');
            fprintf(fout,'\t\t          zout => zpast);\n');
            fprintf(fout,'\n');
        end
        
        fprintf(fout,'\tproc_counter: process(clk,reset)\n');
        fprintf(fout,'\tbegin\n');
        fprintf(fout,'\t\tif reset = ''0'' then\n');
        fprintf(fout,'\t\t\tcounter <= 0;\n');
        fprintf(fout,'\t\telsif rising_edge(clk) then\n');
        fprintf(fout,'\t\t\tif counter = sampling_latency_VS-1 then\n');
        fprintf(fout,'\t\t\t\tcounter <= 0;\n');
        fprintf(fout,'\t\t\telse\n');
        fprintf(fout,'\t\t\t\tcounter <= counter+1;\n');
        fprintf(fout,'\t\t\tend if;\n');
        fprintf(fout,'\t\tend if;\n');
        fprintf(fout,'\tend process;\n');
        fprintf(fout,'\n');
        fprintf(fout,'\tproc_counter_full: process(clk,reset,signal_sample)\n');
        fprintf(fout,'\tbegin\n');
        fprintf(fout,'\t\tif reset = ''0'' then\n');
        fprintf(fout,'\t\t\tcounter_full <= 0;\n');
        fprintf(fout,'\t\telsif rising_edge(clk) and signal_sample = ''1'' then\n');
        fprintf(fout,'\t\t\tif counter_full = MMAX_VS then\n');
        fprintf(fout,'\t\t\t\tcounter_full <= MMAX_VS;\n');
        fprintf(fout,'\t\t\telse\n');
        fprintf(fout,'\t\t\t\tcounter_full <= counter_full+1;\n');
        fprintf(fout,'\t\t\tend if;\n');
        fprintf(fout,'\t\tend if;\n');
        fprintf(fout,'\tend process;\n');
        fprintf(fout,'\n');
        fprintf(fout,'\tsignal_full <= ''1'' when counter_full = MMAX_VS else \n');
        fprintf(fout,'\t                 ''0'';\n');
        fprintf(fout,'\n');
        fprintf(fout,'\tsignal_sample <= ''1'' when counter = 0 else \n');
        fprintf(fout,'\t                 ''0'';\n');
        fprintf(fout,'\n');
        fprintf(fout,'\tsignal_start <= ''1'' when counter = 1 and signal_full = ''1'' else\n');
        fprintf(fout,'\t                ''0'';\n');
        fprintf(fout,'\n');
        % TO DO
        if object.isReducedComplexity
            fprintf(fout,'\tproc_finish: process(clk,reset,');
            for i = 1:npwas
                fprintf(fout,'signal_done%d,',i);
            end
            fprintf(fout,'signal_start)\n');
            fprintf(fout,'\tbegin\n');
            fprintf(fout,'\t\tif reset = ''0'' then\n');
            for i = 1:npwas
                fprintf(fout,'\t\t\tsignal_finish%d <= ''0'';\n',i);
            end
            fprintf(fout,'\t\telsif rising_edge(clk) then\n');
            for i = 1:npwas
                fprintf(fout,'\t\t\tif signal_done%d = ''1'' then\n',i);
                fprintf(fout,'\t\t\t\tsignal_finish%d <= ''1'';\n',i);
                fprintf(fout,'\t\t\tend if;\n');
                fprintf(fout,'\t\t\tif signal_start = ''1'' then\n');
                fprintf(fout,'\t\t\t\tsignal_finish%d <= ''0'';\n',i);
                fprintf(fout,'\t\t\tend if;\n');
            end
            fprintf(fout,'\t\t\tif signal_done = ''1'' then\n');
            for i = 1:npwas
                fprintf(fout,'\t\t\t\tsignal_finish%d <= ''0'';\n',i);
            end
            fprintf(fout,'\t\t\tend if;\n');
            fprintf(fout,'\t\tend if;\n');
            fprintf(fout,'\tend process;\n');
            fprintf(fout,'\t\n');
            fprintf(fout,'\tproc_done: process(clk,reset,signal_finish,signal_done)\n');
            fprintf(fout,'\tbegin\n');
            fprintf(fout,'\t\tif reset = ''0'' then\n');
            fprintf(fout,'\t\t\tsignal_done <= ''0'';\n');
            fprintf(fout,'\t\telsif rising_edge(clk) then\n');
            fprintf(fout,'\t\t\tif signal_finish = ''1'' then\n');
            fprintf(fout,'\t\t\t\tsignal_done <= ''1'';\n');
            fprintf(fout,'\t\t\tend if;\n');
            fprintf(fout,'\t\t\tif signal_done = ''1'' then\n');
            fprintf(fout,'\t\t\t\tsignal_done <= ''0'';\n');
            fprintf(fout,'\t\t\tend if;\n');
            fprintf(fout,'\t\tend if;\n');
            fprintf(fout,'\tend process;\n');
            fprintf(fout,'\t           \n');
            fprintf(fout,'\tsignal_finish <= ');
           for i = 1:npwas-1
               fprintf(fout,'signal_finish%d and ',i);
           end
           fprintf(fout,'signal_finish%d;\n',npwas);
            if strcmp(opts.outputRepresentation,'signed')
               fprintf(fout,'\tsignal_z <= (others => ''0'') when reset = ''0'' else\n');
               fprintf(fout,'\t            std_logic_vector(');
               for i = 1:npwas-1
                   fprintf(fout,'signed(signal_z%d)+',i);
               end
               fprintf(fout,'signed(signal_z%d)) when signal_finish = ''1'' and rising_edge(clk);\n',npwas);
            else
                fprintf(fout,'\tsignal_z <= (others => ''0'') when reset = ''0'' else\n');
                fprintf(fout,'\t            std_logic_vector(');
                for i = 1:npwas-1
                    fprintf(fout,'unsigned(signal_z%d)+',i);
                end
                fprintf(fout,'unsigned(signal_z%d)) when signal_finish = ''1'' and rising_edge(clk);\n',npwas);
            end
           fprintf(fout,'\n');

        else
            fprintf(fout,'\tsignal_done <= signal_done1;\n');
            fprintf(fout,'\tsignal_z <= signal_z1;\n');
        end
        fprintf(fout,'\tz <= signal_z;\n');
        fprintf(fout,'\tdone <= signal_done;\n');
        fprintf(fout,'\n');
        
        fprintf(fout,'\n');
        fprintf(fout,'\tend Behavioral;\n');
        
        
        
        fprintf(fout,'\n');
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'virtualSensor.vhd']);
