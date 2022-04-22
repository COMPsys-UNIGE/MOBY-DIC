function ct_simulinkClosedLoopSim(object,options)
% ct_simulinkClosedLoopSim    Generates a Simulink model for the
%                             simulation of a continuous-time closed-loop
%                            dynamical system
%
% This is a private method called by method ltiSys/generateSimulinkModel.
%
% See also: ltiSys/generateSimulinkModel.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
%
% Copyright (C) 2015 University of Genoa, Italy.

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


simulateVHDL = options.simulateVHDL;
circuit_parameters = options.circuit_parameters;
% Default value
choice = 'Overwrite';

% If the model already exists, ask what to do...
if exist([options.folder,'pwaSys_model.mdl'],'file')
    
    % Construct a questdlg with three options
    choice = questdlg(['A model already exists in ',options.folder,'. Do you want to overwrite it completely or just update the system matrices?'], ...
        ' ', ...
        'Overwrite','Update','Cancel','Cancel');
    
end

% Extract matrices
[A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy,sim_H,sim_K] = object.getMatrices();

% Number of variables
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;
nxref = numel(object.getController.getTrackingVariable);

nDyn = object.nDyn;

sim_A = cell(nDyn,1);
sim_B = cell(nDyn,1);
sim_C = cell(nDyn,1);
sim_D = cell(nDyn,1);

for i=1:nDyn
    % Put together matrices
    sim_A{i} = A{i}; %#ok<*NASGU>
    sim_B{i} = [B{i} Ex{i} Fx{i} Gx{i}];
    sim_C{i} = [eye(object.nx);C{i}];
    sim_D{i} = [zeros(object.nx,object.nu+object.np+object.nd+1);[D{i} Ey{i} Fy{i} Gy{i}]];
end

sim_xnames = object.getStateNames();
sim_unames = object.getInputNames();
sim_pnames = object.getParameterNames();
sim_ynames = object.getOutputNames();
sim_dnames = object.getUnmeasurableInputNames();

% Extract controller
ctrl = object.controller;

if strcmpi(choice,'Overwrite')
    
    preFolder = '';
    % Create folder
    mkdir(options.folder)
    if ~simulateVHDL
    
    % Choose model to generate
    if np > 0
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_model_pd.mdl';
            else
                model = 'pwaSys_reg_model_pd.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_model_p.mdl';
            else
                model = 'pwaSys_reg_model_p.mdl';
            end
        end
    else
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_model_d.mdl';
            else
                model = 'pwaSys_reg_model_d.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_model.mdl';
            else
                model = 'pwaSys_reg_model.mdl';
            end
        end
    end
    else
        % Choose model to generate
    if np > 0
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_model_pd.mdl';
            else
                model = 'pwaSys_reg_model_pd.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_model_p.mdl';
            else
                model = 'pwaSys_reg_model_p.mdl';
            end
        end
    else
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_model_d.mdl';
            else
                model = 'pwaSys_reg_model_d.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_model.mdl';
            else
                model = 'pwaSys_reg_model.mdl';
            end
        end
    end
    end
    
    
    % Copy S-function
    copyfile([options.simulinkFolder,'ct/','s_pwaSys_c.m'],[options.folder,'s_pwaSys_c.m']);
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_np = np;
    sim_nu = nu;
    sim_nxref = nxref;
    sim_ctrl = ctrl;
    
    if ~simulateVHDL
        copyfile([options.simulinkFolder,'ct/s_MPCctrl.m'],[options.folder,'s_MPCctrl.m']);
    
    % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
        'sim_H','sim_K','sim_nx','sim_ny','sim_nd','sim_ctrl',...
        'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
    else
        synthInfo = ctrl.generateVHDL(circuit_parameters);
        
        range = synthInfo.range;
        frequency = synthInfo.circuit_parameters.frequency;
        
        inputRange = synthInfo.circuit_parameters.inputRange;
        
         sim_x_scale_gain = (synthInfo.circuit_parameters.inputRange.max(1:nx)-synthInfo.circuit_parameters.inputRange.min(1:nx))./...
            (range.xmax(:)-range.xmin(:));
        sim_x_scale_bias = (range.xmax(:)'+range.xmin(:)')/2-...
            ((synthInfo.circuit_parameters.inputRange.max(1:nx)+synthInfo.circuit_parameters.inputRange.min(1:nx))/2./sim_x_scale_gain)';
        
        sim_p_scale_gain = (synthInfo.circuit_parameters.inputRange.max(nx+1:nx+np)-synthInfo.circuit_parameters.inputRange.min(nx+1:nx+np))./...
            (range.pmax(:)-range.pmin(:));
        sim_p_scale_bias = (range.pmax(:)'+range.pmin(:)')/2-...
            ((synthInfo.circuit_parameters.inputRange.max(nx+1:nx+np)+synthInfo.circuit_parameters.inputRange.min(nx+1:nx+np))/2./sim_p_scale_gain)';
        
        
        sim_d_scale_gain = (synthInfo.circuit_parameters.inputRange.max(nx+np+1:nx+np+nd)-synthInfo.circuit_parameters.inputRange.min(nx+np+1:nx+np+nd))./...
            (range.dmax(:)-range.dmin(:));
        sim_d_scale_bias = (range.dmax(:)'+range.dmin(:)')/2-...
            ((synthInfo.circuit_parameters.inputRange.max(nx+np+1:nx+np+nd)+synthInfo.circuit_parameters.inputRange.min(nx+np+1:nx+np+nd))/2./sim_d_scale_gain)';
        
        
        sim_xref_scale_gain = (synthInfo.circuit_parameters.inputRange.max(nx+np+nd+1:nx+np+nd+nxref)-synthInfo.circuit_parameters.inputRange.min(nx+np+nd+1:nx+np+nd+nxref))./...
            (range.xrefmax(:)-range.xrefmin(:));
        sim_xref_scale_bias = (range.xrefmax(:)'+range.xrefmin(:)')/2-...
            ((synthInfo.circuit_parameters.inputRange.max(nx+np+nd+1:nx+np+nd+nxref)+synthInfo.circuit_parameters.inputRange.min(nx+np+nd+1:nx+np+nd+nxref))/2./sim_xref_scale_gain)';
        
        
        sim_u_scale_gain = (range.umax(:)-range.umin(:))./(synthInfo.circuit_parameters.outputRange.max-synthInfo.circuit_parameters.outputRange.min);
        sim_u_scale_bias = (range.umax+range.umin)/2-...
            ((synthInfo.circuit_parameters.outputRange.max+synthInfo.circuit_parameters.outputRange.min)/2.*sim_u_scale_gain)';
        
        
        clk_period = 1/frequency;
        sim_n_bit_in = synthInfo.circuit_parameters.inputResolution;
        sim_n_bit_out = synthInfo.circuit_parameters.outputResolution;
        
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D','sim_H','sim_K',...
            'sim_nx','sim_ny','sim_np','sim_nd','sim_nxref','sim_u_scale_gain','sim_u_scale_bias',...
            'sim_x_scale_gain','sim_x_scale_bias','sim_p_scale_gain','sim_p_scale_bias',...
            'sim_d_scale_gain','sim_d_scale_bias','sim_xref_scale_gain','sim_xref_scale_bias',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames','clk_period',...
            'sim_n_bit_in','sim_n_bit_out','sim_nu','sim_nxref');
        
        preFolder = ['circuit_model/',synthInfo.circuit_parameters.inputRepresentation,'_in/'];
    end
    
    % Copy model templates
    copyfile([options.simulinkFolder,'ct/',preFolder,model],[options.folder,'pwaSys_model.mdl']);
    copyfile([options.simulinkFolder,'ct/','plotScope.m'],[options.folder,'plotScope.m']);
    
    
    
    if simulateVHDL
        % Generate black box configuration function
        VHDL_folder = [options.folder,'VHDLfile/'];
        outputRepresentation = synthInfo.circuit_parameters.outputRepresentation;
        
        fp = fopen([options.folder,'embeddedSystem_config.m'],'w');
        
        fprintf(fp,'function embeddedSystem_config(this_block)\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.setTopLevelLanguage(''VHDL'');\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.setEntityName(''switchedControllerInterface'');\n');
        fprintf(fp,'\n');
        fprintf(fp,'\t%% System Generator has to assume that your entity  has a combinational feed through;\n');
        fprintf(fp,'\t%%   if it  doesn''t, then comment out the following line:\n');
        fprintf(fp,'\tthis_block.tagAsCombinational;\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.addSimulinkInport(''reset'');\n');
        for i = 1:nx
            fprintf(fp,'\tthis_block.addSimulinkInport(''x%d'');\n',i);
        end
        for i = 1:np
            fprintf(fp,'\tthis_block.addSimulinkInport(''p%d'');\n',i);
        end
        for i = 1:nd
            fprintf(fp,'\tthis_block.addSimulinkInport(''d%d'');\n',i);
        end
        for i = 1:nxref
            fprintf(fp,'\tthis_block.addSimulinkInport(''xref%d'');\n',i);
        end
        fprintf(fp,'\n');
        for i = 1:nu
            fprintf(fp,'\tthis_block.addSimulinkOutport(''u%d'');\n',i);
        end
        fprintf(fp,'\tthis_block.addSimulinkOutport(''sample'');\n');
        fprintf(fp,'\tthis_block.addSimulinkOutport(''done'');\n');
        
        fprintf(fp,'\n');
        for i = 1:nu
            fprintf(fp,'\tu%d_port = this_block.port(''u%d'');\n',i,i);
            if strcmp(outputRepresentation,'signed')
                fprintf(fp,'\tu%d_port.setType(''Fix_%d_0'');\n',i,sim_n_bit_out);
            else
                fprintf(fp,'\tu%d_port.setType(''uFix_%d_0'');\n',i,sim_n_bit_out);
            end
        end
        fprintf(fp,'\tsample_port = this_block.port(''sample'');\n');
        fprintf(fp,'\tsample_port.setType(''UFix_1_0'');\n');
        fprintf(fp,'\tsample_port.useHDLVector(false);\n');
        fprintf(fp,'\tdone_port = this_block.port(''done'');\n');
        fprintf(fp,'\tdone_port.setType(''UFix_1_0'');\n');
        fprintf(fp,'\tdone_port.useHDLVector(false);\n');
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
        for i = 1:nx
            fprintf(fp,'\t\tif (this_block.port(''x%d'').width ~= %d);\n',i,sim_n_bit_in);
            fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "x%d" must have width=%d.'');\n',i,sim_n_bit_in);
            fprintf(fp,'\t\tend\n');
            fprintf(fp,'\n');
        end
        for i = 1:np
            fprintf(fp,'\t\tif (this_block.port(''p%d'').width ~= %d);\n',i,sim_n_bit_in);
            fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "p%d" must have width=%d.'');\n',i,sim_n_bit_in);
            fprintf(fp,'\t\tend\n');
            fprintf(fp,'\n');
        end
        for i = 1:nd
            fprintf(fp,'\t\tif (this_block.port(''d%d'').width ~= %d);\n',i,sim_n_bit_in);
            fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "d%d" must have width=%d.'');\n',i,sim_n_bit_in);
            fprintf(fp,'\t\tend\n');
            fprintf(fp,'\n');
        end
        for i = 1:nxref
            fprintf(fp,'\t\tif (this_block.port(''xref%d'').width ~= %d);\n',i,sim_n_bit_in);
            fprintf(fp,'\t\t\tthis_block.setError(''Input data type for port "xref%d" must have width=%d.'');\n',i,sim_n_bit_in);
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
        %find all VHDL files
        MyDirInfo = dir([VHDL_folder,'*.vhd*']);
        for i=1:numel(MyDirInfo);
            fprintf(fp,'this_block.addFile(''%s'');\n',['VHDLfile/',MyDirInfo(i).name]);
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
    
    % Generate configuration file
    fp = fopen([options.folder,'setParameters.m'],'w');
    
    fprintf(fp,'%% SET PARAMETERS\n');
    fprintf(fp,'%%\n');
    fprintf(fp,'%% In this script you can set the parameters for the Simulink simulation of\n');
    fprintf(fp,'%% the open-loop continuous-time system\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% Initial condition\n');
    fprintf(fp,'sim_x0 = [');
    for i = 1:object.nx-1
        fprintf(fp,'0;');
    end
    fprintf(fp,'0];\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% Simulation time\n');
    fprintf(fp,'sim_T = 1;\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% If you do not use the signal builder for the inputs, parameters and \n');
    fprintf(fp,'%% unmeasurable inputs, you can set their constant values below\n');
    fprintf(fp,'\n');
    
    if np > 0
        fprintf(fp,'%% Constant parameters\n');
        fprintf(fp,'sim_p = [');
        for i = 1:object.np-1
            fprintf(fp,'0;');
        end
        fprintf(fp,'0];\n');
        fprintf(fp,'\n');
    end
    if nd > 0
        fprintf(fp,'%% Constant unmeasurable inputs\n');
        fprintf(fp,'sim_d = [');
        for i = 1:object.nd-1
            fprintf(fp,'0;');
        end
        fprintf(fp,'0];\n');
        fprintf(fp,'\n');
    end
    if ctrl.isTracking
        fprintf(fp,'%% Reference state\n');
        fprintf(fp,'sim_xref = 0;\n');
    end
    fclose(fp);
    
    disp('Model succesfully created in:')
    disp(options.folder)
    disp(' ')
    disp('The model will now be opened...')
    disp(' ')
    cd(options.folder)
    edit('setParameters.m')
    pwaSys_model
    
    helpdlg('Before running the model, set the parameters in file "setParameters.m"')
    
elseif strcmpi(choice,'Update')
    
     if ~simulateVHDL
        error('VHDL co-simulation model cannor be updated');
    end
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_ctrl = ctrl;
    
    try
        % Save .mat files with the matrices
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
            'sim_H','sim_K','sim_nx','sim_ny','sim_nd','sim_ctrl',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
        
    catch
        error('Model cannot be updated. Maybe you are in the wrong folder.');
    end
    
    disp('Model succesfully updated.')
    disp(' ')
    disp('The model will now be opened...')
    disp(' ')
    cd(options.folder)
    ltiSys_model
    
else
    
    disp('Canceled.')
    disp(' ')
    
end






