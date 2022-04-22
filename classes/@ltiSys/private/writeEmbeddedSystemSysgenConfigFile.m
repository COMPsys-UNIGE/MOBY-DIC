function writeEmbeddedSystemSysgenConfigFile(object, options)
% writeEmbeddedSystemSysgenConfigFile    
% Writes the System Generator configuration file for the embedded system
% (LTI system + observer + controller)
%
% This is a private method called by method ltiSys/generateSimulinkModel.
%
% writeEmbeddedSystemSysgenConfigFile(OBJ,OPTS)
% OBJ is a ltiSys object, OPTS is the struct containing the options for the
% generation of the VHDL files of the controller associated to OBJ.
%
% See also: ltiSys/generateSimulinkModel.

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

nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;

obs = object.getObserver;
ctrl = object.getController;

nxref = numel(ctrl.getTrackingVariable);  
        
range = options.circuit_parameters.range;
ctrlPar = MPCctrlVHDLset(ctrl, options.circuit_parameters);
ctrlPar.generateInterface = false;
                
obsPar = observerVHDLset(obs, options.circuit_parameters);
obsPar.range = range;
obsPar.outputResolution = ctrlPar.inputResolution;
obsPar.outputRange = ctrlPar.inputRange;
obsPar.outputRepresentation = ctrlPar.inputRepresentation;
obsPar.generateInterface = false;

ctrlPar.inputResolution = obsPar.outputResolution;
ctrlPar.inputRange = obsPar.outputRange;
ctrlPar.inputRepresentation = obsPar.inputRepresentation;

VHDL_folder = options.circuit_parameters.folder;

outputRepresentation = ctrlPar.outputRepresentation;

sim_n_bit_in = obsPar.inputResolution;
sim_n_bit_out = ctrlPar.outputResolution;

fp = fopen([options.folder,'sysgen_config.m'],'w');

fprintf(fp,'function sysgen_config(this_block)\n');
fprintf(fp,'\n');
fprintf(fp,'\tthis_block.setTopLevelLanguage(''VHDL'');\n');
fprintf(fp,'\n');
if isa(ctrl, 'implicitMPCctrl')
    fprintf(fp,'\tthis_block.setEntityName(''embeddedSystemInterface'');\n');
else
    fprintf(fp,'\tthis_block.setEntityName(''embeddedSystem'');\n');
end
fprintf(fp,'\n');
fprintf(fp,'\t%% System Generator has to assume that your entity  has a combinational feed through;\n');
fprintf(fp,'\t%%   if it  doesn''t, then comment out the following line:\n');
fprintf(fp,'\tthis_block.tagAsCombinational;\n');
fprintf(fp,'\n');

fprintf(fp,'\tthis_block.addSimulinkInport(''reset'');\n');

if isa(ctrl, 'implicitMPCctrl')
    fprintf(fp,'\tthis_block.addSimulinkInport(''obs_start'');\n');
    fprintf(fp,'\tthis_block.addSimulinkInport(''ctrl_start'');\n');
end
for i = 1:np
    fprintf(fp,'\tthis_block.addSimulinkInport(''p%d'');\n',i);
end
for i = 1:nxref
    fprintf(fp,'\tthis_block.addSimulinkInport(''xref%d'');\n',i);
end
for i = 1:ny
    fprintf(fp,'\tthis_block.addSimulinkInport(''y%d'');\n',i);
end
fprintf(fp,'\n');
for i = 1:nu
    fprintf(fp,'\tthis_block.addSimulinkOutport(''u%d'');\n',i);
end
fprintf(fp,'\tthis_block.addSimulinkOutport(''sample'');\n');
fprintf(fp,'\tthis_block.addSimulinkOutport(''output_ready'');\n');

fprintf(fp,'\n');
for i = 1:nu
    fprintf(fp,'\tu%d_port = this_block.port(''u%d'');\n',i,i);
    
        if isa(ctrl,'implicitMPCctrl')
            if ctrlPar.useDAC==1
                sim_bin_point_out = 0;
            else
                sim_bin_point_out = sim_n_bit_out - ceil(log2(max(abs([options.range.umax, options.range.umin]),[],'all')));
            end
        else
            sim_bin_point_out = 0;
        end
        
        if strcmp(outputRepresentation,'signed')
            fprintf(fp,'\tu%d_port.setType(''Fix_%d_%d'');\n', i, sim_n_bit_out, max(0,sim_bin_point_out-1));
        else
            fprintf(fp,'\tu%d_port.setType(''UFix_%d_%d'');\n', i, sim_n_bit_out, sim_bin_point_out);
        end
end

fprintf(fp,'\tsample_port = this_block.port(''sample'');\n');
fprintf(fp,'\tsample_port.setType(''UFix_1_0'');\n');
fprintf(fp,'\tsample_port.useHDLVector(false);\n');
fprintf(fp,'\toutput_ready_port = this_block.port(''output_ready'');\n');
fprintf(fp,'\toutput_ready_port.setType(''UFix_1_0'');\n');
fprintf(fp,'\toutput_ready_port.useHDLVector(false);\n');
fprintf(fp,'\n');
fprintf(fp,'\t%% -----------------------------\n');
fprintf(fp,'\tif (this_block.inputTypesKnown)\n');
fprintf(fp,'\t\t%% do input type checking, dynamic output type and generic setup in this code block.\n');
fprintf(fp,'\n');
fprintf(fp,'\t\tif (this_block.port(''reset'').width ~= 1);\n');
fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "reset" must have width=1.'');\n');
fprintf(fp,'\t\tend\n');
fprintf(fp,'\n');
fprintf(fp,'\t\tthis_block.port(''reset'').useHDLVector(false);\n');
fprintf(fp,'\n');

if isa(ctrl, 'implicitMPCctrl')
    fprintf(fp,'\t\tif (this_block.port(''obs_start'').width ~= 1);\n');
    fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "obs_start" must have width=1.'');\n');
    fprintf(fp,'\t\tend\n');
    fprintf(fp,'\n');
    fprintf(fp,'\t\tthis_block.port(''obs_start'').useHDLVector(false);\n');
    fprintf(fp,'\n');
    
    fprintf(fp,'\t\tif (this_block.port(''ctrl_start'').width ~= 1);\n');
    fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "ctrl_start" must have width=1.'');\n');
    fprintf(fp,'\t\tend\n');
    fprintf(fp,'\n');
    fprintf(fp,'\t\tthis_block.port(''ctrl_start'').useHDLVector(false);\n');
    fprintf(fp,'\n');
end

for i = 1:np
    fprintf(fp,'\t\tif (this_block.port(''p%d'').width ~= %d);\n',i,sim_n_bit_in);
    fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "p%d" must have width=%d.'');\n',i,sim_n_bit_in);
    fprintf(fp,'\t\tend\n');
    fprintf(fp,'\n');
end

for i = 1:nxref
    fprintf(fp,'\t\tif (this_block.port(''xref%d'').width ~= %d);\n',i,sim_n_bit_in);
    fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "xref%d" must have width=%d.'');\n',i,sim_n_bit_in);
    fprintf(fp,'\t\tend\n');
    fprintf(fp,'\n');
end

for i = 1:ny
    fprintf(fp,'\t\tif (this_block.port(''y%d'').width ~= %d);\n',i,sim_n_bit_in);
    fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "y%d" must have width=%d.'');\n',i,sim_n_bit_in);
    fprintf(fp,'\t\tend\n');
    fprintf(fp,'\n');
end

fprintf(fp,'\tend  %% if(inputTypesKnown)\n');
fprintf(fp,'\t%% -----------------------------\n');
fprintf(fp,'\n');
fprintf(fp,'\t%% -----------------------------\n');
fprintf(fp,'\tif (this_block.inputRatesKnown)\n');
fprintf(fp,'\t\tsetup_as_single_rate(this_block,''clk'',''ce'')\n');
fprintf(fp,'\tend  %% if(inputRatesKnown)\n');
fprintf(fp,'\t%% -----------------------------\n');
fprintf(fp,'\n');
fprintf(fp,'\t%% (!) Set the inout port rate to be the same as the first input\n');
fprintf(fp,'\t%%     rate. Change the following code if this is untrue.\n');
fprintf(fp,'\tuniqueInputRates = unique(this_block.getInputRates);\n');
fprintf(fp,'\n');
fprintf(fp,'\n');
fprintf(fp,'\t%% Add addtional source files as needed.\n');
fprintf(fp,'\t%%  |-------------\n');
fprintf(fp,'\t%%  | Add files in the order in which they should be compiled.\n');
fprintf(fp,'\t%%  | If two files "a.vhd" and "b.vhd" contain the entities\n');
fprintf(fp,'\t%%  | entity_a and entity_b, and entity_a contains a\n');
fprintf(fp,'\t%%  | component of type entity_b, the correct sequence of\n');
fprintf(fp,'\t%%  | addFile() calls would be:\n');
fprintf(fp,'\t%%  |    this_block.addFile(''b.vhd'');\n');
fprintf(fp,'\t%%  |    this_block.addFile(''a.vhd'');\n');
fprintf(fp,'\t%%  |-------------\n');
fprintf(fp,'\n');
fprintf(fp,'\n');

% Find all VHDL files
if ~isa(ctrl, 'implicitMPCctrl')
    MyDirInfo = dir([VHDL_folder,'/*.vhd*']);
    for i=1:numel(MyDirInfo)
        fprintf(fp,'this_block.addFile(''%s'');\n',['VHDLfile/',MyDirInfo(i).name]);
    end
else
    MyDirInfo = dir([VHDL_folder,'/vhdl/','*.vhd*']);
    for i=1:numel(MyDirInfo)
        fprintf(fp,'\tthis_block.addFile(''%s'');\n',['VHDLfile/vhdl/',MyDirInfo(i).name]);
    end
end


fprintf(fp,'\n');
fprintf(fp,'return;\n');
fprintf(fp,'\n');
fprintf(fp,'\n');
fprintf(fp,'%% ------------------------------------------------------------\n');
fprintf(fp,'\n');
fprintf(fp,'function setup_as_single_rate(block,clkname,cename)\n');
fprintf(fp,'\tinputRates = block.inputRates;\n');
fprintf(fp,'\tuniqueInputRates = unique(inputRates);\n');
fprintf(fp,'\tif (length(uniqueInputRates)==1 & uniqueInputRates(1)==Inf)\n');
fprintf(fp,'\t\tblock.addError(''The inputs to this block cannot all be constant.'');\n');
fprintf(fp,'\t\treturn;\n');
fprintf(fp,'\tend\n');
fprintf(fp,'\tif (uniqueInputRates(end) == Inf)\n');
fprintf(fp,'\t\thasConstantInput = true;\n');
fprintf(fp,'\t\tuniqueInputRates = uniqueInputRates(1:end-1);\n');
fprintf(fp,'\tend\n');
fprintf(fp,'\tif (length(uniqueInputRates) ~= 1)\n');
fprintf(fp,'\t\tblock.addError(''The inputs to this block must run at a single rate.'');\n');
fprintf(fp,'\t\treturn;\n');
fprintf(fp,'\tend\n');
fprintf(fp,'\ttheInputRate = uniqueInputRates(1);\n');
fprintf(fp,'\tfor i = 1:block.numSimulinkOutports\n');
fprintf(fp,'\t\tblock.outport(i).setRate(theInputRate);\n');
fprintf(fp,'\tend\n');
fprintf(fp,'\tblock.addClkCEPair(clkname,cename,theInputRate);\n');
fprintf(fp,'\treturn;\n');
fprintf(fp,'\n');
fprintf(fp,'%% ------------------------------------------------------------\n');
fclose(fp);

end

