%% WRITE TEST
% Writes the VHDL code implementing block test_vs_ser_interface or block
% test_vs_par_interface in vs circuit
%
% SYNTAX
%
% writeTest(object,synthesisInfo,folder)
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

function writeTest(object,synthesisInfo,folder)

type = synthesisInfo.type;

nu = object.getNumberOfInputs();
ny = object.getNumberOfMeasurableOutputs();
mu = object.getInputTimeWindow();
my = object.getInputTimeWindow();

current = object.isCurrent();

if strcmpi(type,'serial')
    fout = fopen([folder,'test_vs_ser_interface.vhd'],'w');
    entity_name = 'vs_ser_interface';
else
    fout = fopen([folder,'test_vs_par_interface.vhd'],'w');
    entity_name = 'vs_par_interface';
end

fprintf(fout,'LIBRARY IEEE;\n');
fprintf(fout,'USE IEEE.STD_LOGIC_1164.ALL;\n');
fprintf(fout,'USE IEEE.NUMERIC_STD.ALL;\n');
fprintf(fout,'\n');
fprintf(fout,'-- Necessary for reading and writing files\n');
fprintf(fout,'LIBRARY STD;\n');
fprintf(fout,'USE STD.TEXTIO.ALL;\n');
fprintf(fout,'\n');
if strcmpi(type,'serial')
    fprintf(fout,'use work.vs_ser_package.all;\n');
else
    fprintf(fout,'use work.vs_par_package.all;\n');
end
fprintf(fout,'\n');
fprintf(fout,'-- Entity declaration\n');
fprintf(fout,'\n');
fprintf(fout,'ENTITY test_%s IS\n',entity_name);
fprintf(fout,'\n');
fprintf(fout,'END test_%s;\n',entity_name);
fprintf(fout,'\n');
fprintf(fout,'-- Architecture declaration\n');
fprintf(fout,' \n');
fprintf(fout,'ARCHITECTURE behavior OF test_%s IS\n',entity_name);
fprintf(fout,'\n');
    fprintf(fout,'-- Component Declaration for the Unit Under Test (UUT)\n');
fprintf(fout,'\n');
	
fprintf(fout,'COMPONENT %s\n',entity_name);
fprintf(fout,'\tPort ( clk : in std_logic;\n');
fprintf(fout,'\t\tce : in std_logic;\n');
fprintf(fout,'\t\treset : in std_logic;\n');

for j = 1:nu
    fprintf(fout,'\t\tu%d : in  STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n',j);
end

for j = 1:ny
    fprintf(fout,'\t\ty%d : in  STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n',j);
end

fprintf(fout,'\t\tz : out  STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n');
fprintf(fout,'\t\tready : out  STD_LOGIC);\n');
fprintf(fout,'END COMPONENT;\n');
fprintf(fout,'\n');
fprintf(fout,'-- Signals declaration\n');
fprintf(fout,'\n');

for j = 1:nu
    fprintf(fout,'\tsignal u%d : STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n',j);
end

for j = 1:ny
    fprintf(fout,'\tsignal y%d : STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n',j);
end

fprintf(fout,'\tsignal z : STD_LOGIC_VECTOR (N_BIT-1 downto 0);\n');

fprintf(fout,'\tsignal clk : std_logic := ''0'';\n');
fprintf(fout,'\tsignal reset : std_logic := ''0'';\n');
fprintf(fout,'\tsignal ready : std_logic;\n');
fprintf(fout,'\n');
fprintf(fout,'\t-- Delay of input signals\n');
fprintf(fout,'\tconstant delay : time := CLK_PERIOD/4;\n');
fprintf(fout,'\n');
fprintf(fout,'-- Sampling time of the circuit\n');
fprintf(fout,'constant SAMPLING_TIME: time := (SAMPLING_LATENCY-ACQUISITION_LATENCY-1)*CLK_PERIOD;\n');
fprintf(fout,'\n');
fprintf(fout,'-- When ''1'' the pwas function has been written to file, preventing\n');
fprintf(fout,'-- multiple writings\n');
fprintf(fout,'signal write_performed : std_logic;\n');
fprintf(fout,'\n');
fprintf(fout,'--Files\n');
fprintf(fout,'FILE input_file :text open READ_MODE is "./test/inputs.txt";\n');
fprintf(fout,'FILE output_file :text open WRITE_MODE is "./test/outputs.txt"; \n');
fprintf(fout,'\n');

fprintf(fout,'BEGIN\n');

fprintf(fout,'\t-- Instantiate the Unit Under Test (UUT)\n');
fprintf(fout,'\tuut: %s\n',entity_name);
fprintf(fout,'\tPORT MAP (\n');
fprintf(fout,'\t\tclk => clk,\n');
fprintf(fout,'\t\tce => ''1'',\n');
fprintf(fout,'\t\treset => reset,\n');
for j = 1:nu
    fprintf(fout,'\t\tu%d => u%d,\n',j,j);
end
for j = 1:ny
    fprintf(fout,'\t\ty%d => y%d,\n',j,j);
end
fprintf(fout,'\t\tz => z,\n');
fprintf(fout,'\t\tready => ready\n');
fprintf(fout,'\t);\n');

fprintf(fout,'\n');

fprintf(fout,'\t-- Clock process\n');
fprintf(fout,'\tproc_clock :process\n');
fprintf(fout,'\tbegin\n');
fprintf(fout,'\tclk <= ''0'';\n');
fprintf(fout,'\twait for CLK_PERIOD/2;\n');
fprintf(fout,'\tclk <= ''1'';\n');
fprintf(fout,'\twait for CLK_PERIOD/2;\n');
fprintf(fout,'\tend process;\n');
fprintf(fout,'\t\n');
fprintf(fout,'\t-- Reset process\n');
fprintf(fout,'\tproc_reset :process\n');
fprintf(fout,'\tbegin\n');
fprintf(fout,'\t-- Activate reset for 3 clock cycles\n');
fprintf(fout,'\treset <= ''0'';\n');
fprintf(fout,'\twait for 3*CLK_PERIOD;\n');
fprintf(fout,'\t\n');
fprintf(fout,'\t-- Turn off reset\n');
fprintf(fout,'\treset <= ''1'';\n');
fprintf(fout,'\t\n');
fprintf(fout,'\twait;\n');
fprintf(fout,'\tend process;\n');

fprintf(fout,'\n');

fprintf(fout,'\t-- Stimulus process\n');
fprintf(fout,'\tproc_read_input: PROCESS\n');

fprintf(fout,'\n');

fprintf(fout,'\tVARIABLE lnin: line;\n');
fprintf(fout,'\tVARIABLE input_x: bit_vector(N_BIT-1 downto 0);\n');
fprintf(fout,'\t\n');
fprintf(fout,'\tBEGIN\n');
fprintf(fout,'\t\n');
fprintf(fout,'\tif reset = ''0'' then\n');

for j = 1:nu
    fprintf(fout,'\t\tu%d <= (others => ''0'');\n',j);
end
for j = 1:ny
    fprintf(fout,'\t\ty%d <= (others => ''0'');\n',j);
end

fprintf(fout,'\t\twait until reset = ''1'';\n');
fprintf(fout,'\t\tend if;\n');

fprintf(fout,'\n');

fprintf(fout,'\t\t-- If reset = ''1'', read input from file\n');
fprintf(fout,'\t\tif not endfile(input_file) then\n');
fprintf(fout,'\t\twait until rising_edge(clk);\n');
fprintf(fout,'\t\twait for delay;\n');

for j = 1:nu
    fprintf(fout,'\t\treadline(input_file,lnin);\n');
    fprintf(fout,'\t\tread(lnin,input_x);\n');
    fprintf(fout,'\t\tu%d <= to_stdlogicvector(input_x);\n',j);
end
            
for j = 1:ny
    fprintf(fout,'\t\treadline(input_file,lnin);\n');
    fprintf(fout,'\t\tread(lnin,input_x);\n');
    fprintf(fout,'\t\ty%d <= to_stdlogicvector(input_x);\n',j);
end

fprintf(fout,'\telse\n');
fprintf(fout,'\t\twait until rising_edge(clk);\n');
fprintf(fout,'\t\twait for delay;\n');
fprintf(fout,'\t\t-- If the file is finished, throw an error message\n');
fprintf(fout,'\t\tassert false report "NONE. End of simulation." severity failure;\n');
fprintf(fout,'\tend if;\n');

fprintf(fout,'\twait for SAMPLING_TIME;\n');

fprintf(fout,'\tEND PROCESS;\n');
fprintf(fout,'\n');
fprintf(fout,'\t-- When ready is one and there is a rising edge of the clock, the output\n');
fprintf(fout,'\t-- (i.e. the value of pwas function) is written to the output file.\n');
fprintf(fout,'\tproc_write_output: PROCESS(z,ready,clk)\n');
fprintf(fout,'\tVARIABLE lnoutf: line;\n');
fprintf(fout,'\tBEGIN\n');
fprintf(fout,'\tif rising_edge(clk) then\n');
fprintf(fout,'\tif ready = ''1'' and write_performed = ''0'' then \n');
fprintf(fout,'\twrite(lnoutf,to_bitvector(z));\n');
fprintf(fout,'\twriteline(output_file,lnoutf);\n');
fprintf(fout,'\tend if;\n');
fprintf(fout,'\tend if;\n');
fprintf(fout,'\tEND PROCESS;\n');

fprintf(fout,'\t-- When the value of the pwas function has been written to file\n');
fprintf(fout,'\t-- write_performed is set to ''1'', so as to avoid multiple writings\n');
fprintf(fout,'\tproc_write_performed: PROCESS\n');

fprintf(fout,'\tBEGIN\n');

fprintf(fout,'\twrite_performed <= ''0'';\n');

fprintf(fout,'\t		wait until ready = ''1'';\n');
fprintf(fout,'\twait for 2*CLK_PERIOD;\n');

fprintf(fout,'\twrite_performed <= ''1'';\n');

fprintf(fout,'\tif ready = ''0'' then\n');
fprintf(fout,'\twrite_performed <= ''0'';\n');
fprintf(fout,'\telse\n');
fprintf(fout,'\twait until ready = ''0'';\n');
fprintf(fout,'\tend if;\n');

fprintf(fout,'\tEND PROCESS;\n');

fprintf(fout,'\tEND;\n');

fclose(fout);

if strcmpi(type,'serial')
    disp(['Generated file ',folder,'test_vs_ser_interface.vhd']);
else
    disp(['Generated file ',folder,'test_vs_par_interface.vhd']);
end