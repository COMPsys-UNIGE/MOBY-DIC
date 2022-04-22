function dt_simulinkClosedLoopObsSim(object,options)
% dt_simulinkClosedLoopObsSim    Generates a Simulink model for the 
%                                simulation of a discrete-time closed-loop 
%                                dynamical system with observer
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

% Extract observer
obs = object.getObserver();

% Extract observer matrices, gain and sampling time
[Aobs, Bobs, Cobs, Dobs, Gxobs, Gyobs,sim_Hobs,sim_Kobs] = obs.getMatrices();

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

% Put together observer matrices
sim_Aobs =  cell(nDyn,1);
sim_Bobs =  cell(nDyn,1);
sim_Cobs = cell(nDyn,1);
sim_Dobs =  cell(nDyn,1);

sim_Cobs_fake = cell(nDyn,1);
sim_Dobs_fake = cell(nDyn,1);

for i=1:nDyn
    % Put together matrices
    sim_A{i} = A{i}; %#ok<*NASGU>
    sim_B{i} = [B{i} Ex{i} Fx{i} Gx{i}];
    sim_C{i} = [eye(object.nx);C{i}];
    sim_D{i} = [zeros(object.nx,object.nu+object.np+object.nd+1);[D{i} Ey{i} Fy{i} Gy{i}]];

    % Put together observer matrices
sim_Aobs{i} = Aobs{i};
sim_Bobs{i} = [Bobs{i} Gxobs{i}];
sim_Cobs{i} = Cobs{i};
sim_Dobs{i} = [Dobs{i} Gyobs{i}];

%fake obser matrices to avoid loop in simulink
sim_Cobs_fake{i} = eye(object.nx+object.nd);
sim_Dobs_fake{i} = zeros(object.nx+object.nd,object.nu+object.np+object.ny+1);
end

sim_xnames = object.getStateNames();
sim_unames = object.getInputNames();
sim_pnames = object.getParameterNames();
sim_ynames = object.getOutputNames();
sim_dnames = object.getUnmeasurableInputNames();

% Extract controller
ctrl = object.controller;

% Sampling time of the observer
sim_TsObs = obs.getSamplingTime();

% System sampling time
sim_Ts = object.getSamplingTime();

if strcmpi(choice,'Overwrite')
    
    preFolder = '';
    % Create folder
    mkdir(options.folder)
    if ~simulateVHDL
    
    % Choose model to generate
    if np > 0
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_pd.mdl';
            else
                model = 'pwaSys_reg_obs_model_pd.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_p.mdl';
            else
                model = 'pwaSys_reg_obs_model_p.mdl';
            end
        end
    else
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_d.mdl';
            else
                model = 'pwaSys_reg_obs_model_d.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model.mdl';
            else
                model = 'pwaSys_reg_obs_model.mdl';
            end
        end
    end
    else
        % Choose model to generate
    if np > 0
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_pd.mdl';
            else
                model = 'pwaSys_reg_obs_model_pd.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_p.mdl';
            else
                model = 'pwaSys_reg_obs_model_p.mdl';
            end
        end
    else
        if nd > 0
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model_d.mdl';
            else
                model = 'pwaSys_reg_obs_model_d.mdl';
            end
        else
            if ctrl.isTracking
                model = 'pwaSys_track_obs_model.mdl';
            else
                model = 'pwaSys_reg_obs_model.mdl';
            end
        end
    end
    end
    
     sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_np = np;
    sim_nu = nu;
    sim_nxref = nxref;
    sim_ctrl = ctrl;
    
       
     
    if ~simulateVHDL 
          % Copy S-function
    copyfile([options.simulinkFolder,'dt/s_MPCctrl.m'],[options.folder,'s_MPCctrl.m']);
    copyfile([options.simulinkFolder,'dt/','pwaMult.m'],[options.folder,'pwaMult.m']);

        
    % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
        'sim_Aobs','sim_Bobs','sim_Cobs','sim_Dobs','sim_TsObs','sim_Ts',...
        'sim_H','sim_K','sim_Hobs','sim_Kobs','sim_Cobs_fake','sim_Dobs_fake',...
        'sim_nx','sim_nu','sim_np','sim_ny','sim_nd','ctrl',...
        'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
   else
     es = embeddedSystem(object,circuit_parameters.range);
        synthInfo = es.generateVHDL(circuit_parameters);
        range = synthInfo.range;
        
        
        frequency = synthInfo.circuit_parameters.frequency;
        
        inputRange = synthInfo.circuit_parameters.inputRange;
        
        sim_y_scale_gain = (synthInfo.circuit_parameters.inputRange.max(1:ny)-synthInfo.circuit_parameters.inputRange.min(1:ny))./...
            (range.ymax(:)-range.ymin(:));
        sim_y_scale_bias = (range.ymax+range.ymin)/2-...
            ((synthInfo.circuit_parameters.inputRange.max(1:ny)+synthInfo.circuit_parameters.inputRange.min(1:ny))/2./sim_y_scale_gain);
        
        if np == 0
            sim_p_scale_gain = [];
            sim_p_scale_bias = [];
        else
            sim_p_scale_gain = (synthInfo.circuit_parameters.inputRange.max(ny+1:np+ny,:)-synthInfo.circuit_parameters.inputRange.min(ny+1:np+ny,:))./...
                (range.pmax(:)-range.pmin(:));
            sim_p_scale_bias = (range.pmax+range.pmin)/2-...
                ((synthInfo.circuit_parameters.inputRange.max(ny+1:np+ny)+synthInfo.circuit_parameters.inputRange.min(ny+1:np+ny))/2./sim_p_scale_gain)';
        end
        
        if nxref == 0
            sim_xref_scale_gain = [];
            sim_xref_scale_bias = [];
        else
            sim_xref_scale_gain = (synthInfo.circuit_parameters.inputRange.max(np+ny+1:np+ny+nxref)-synthInfo.circuit_parameters.inputRange.min(np+ny+1:np+ny+nxref))./...
                (range.xrefmax(:)-range.xrefmin(:));
            sim_xref_scale_bias = (range.xrefmax+range.xrefmin)/2-...
                ((synthInfo.circuit_parameters.inputRange.max(np+ny+1:np+ny+nxref)+synthInfo.circuit_parameters.inputRange.min(np+ny+1:np+ny+nxref))/2./sim_xref_scale_gain)';
        end
        
        sim_u_scale_gain = (range.umax(:)-range.umin(:))./(synthInfo.circuit_parameters.outputRange.max-synthInfo.circuit_parameters.outputRange.min);
        sim_u_scale_bias = (range.umax+range.umin)/2-...
            ((synthInfo.circuit_parameters.outputRange.max+synthInfo.circuit_parameters.outputRange.min)/2.*sim_u_scale_gain)';
        
        
        clk_period = 1/frequency;
        sim_n_bit_in = synthInfo.circuit_parameters.inputResolution;
        sim_n_bit_out = synthInfo.circuit_parameters.outputResolution;
        
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D','sim_H','sim_K','sim_Ts',...
            'sim_nx','sim_ny','sim_np','sim_nd','sim_nxref','sim_u_scale_gain','sim_u_scale_bias',...
            'sim_y_scale_gain','sim_y_scale_bias','sim_p_scale_gain','sim_p_scale_bias',...
            'sim_xref_scale_gain','sim_xref_scale_bias',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames','clk_period',...
            'sim_n_bit_in','sim_n_bit_out','sim_nu','sim_nxref');
        
        preFolder = ['circuit_model/',preFolder,synthInfo.circuit_parameters.inputRepresentation,'_in/'];
        
    end
    
       % Copy model templates
    copyfile([options.simulinkFolder,'dt/',preFolder,model],[options.folder,'pwaSys_model.mdl']);
    copyfile([options.simulinkFolder,'dt/','plotScope.m'],[options.folder,'plotScope.m']); 
      copyfile([options.simulinkFolder,'dt/','s_pwaSys_d.m'],[options.folder,'s_pwaSys_d.m']);

    
       if simulateVHDL
        % Generate black box configuration function
        VHDL_folder = [options.folder,'VHDLfile/'];
        outputRepresentation = synthInfo.circuit_parameters.outputRepresentation;
        
        fp = fopen([options.folder,'embeddedSystem_config.m'],'w');
        
        fprintf(fp,'function embeddedSystem_config(this_block)\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.setTopLevelLanguage(''VHDL'');\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.setEntityName(''embeddedSystem'');\n');
        fprintf(fp,'\n');
        fprintf(fp,'\t%% System Generator has to assume that your entity  has a combinational feed through;\n');
        fprintf(fp,'\t%%   if it  doesn''t, then comment out the following line:\n');
        fprintf(fp,'\tthis_block.tagAsCombinational;\n');
        fprintf(fp,'\n');
        fprintf(fp,'\tthis_block.addSimulinkInport(''reset'');\n');
        for i = 1:np
            fprintf(fp,'\tthis_block.addSimulinkInport(''p%d'');\n',i);
        end
        for i = 1:nxref
            fprintf(fp,'\tthis_block.addSimulinkInport(''x_ref%d'');\n',i);
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
            if strcmp(outputRepresentation,'signed')
                fprintf(fp,'\tu%d_port.setType(''Fix_%d_0'');\n',i,sim_n_bit_out);
            else
                fprintf(fp,'\tu%d_port.setType(''uFix_%d_0'');\n',i,sim_n_bit_out);
            end
        end
        fprintf(fp,'\tsample_port = this_block.port(''sample'');\n');
        fprintf(fp,'\tsample_port.setType(''UFix_1_0'');\n');
        fprintf(fp,'\tsample_port.useHDLVector(false);\n');
        fprintf(fp,'\tdone_port = this_block.port(''output_ready'');\n');
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
        for i = 1:ny
            fprintf(fp,'\t\tif (this_block.port(''y%d'').width ~= %d);\n',i,sim_n_bit_in);
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
        for i = 1:nxref
            fprintf(fp,'\t\tif (this_block.port(''x_ref%d'').width ~= %d);\n',i,sim_n_bit_in);
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
    fprintf(fp,'%% Observer initial condition for the states\n');
    fprintf(fp,'sim_x0obs = [');
    for i = 1:object.nx-1
        fprintf(fp,'0;');
    end
    fprintf(fp,'0];\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% Observer initial condition for the unmeasurable inputs\n');
    fprintf(fp,'sim_d0obs = [');
    for i = 1:object.nd-1
        fprintf(fp,'0;');
    end
    if object.nd > 0
        fprintf(fp,'0];\n');
    else
        fprintf(fp,'];\n');
    end
    fprintf(fp,'\n');
    fprintf(fp,'%% Simulation time\n');
    fprintf(fp,'sim_T = 1;\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% If you do not use the signal builder for the parameters and \n');
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
    sim_nu = nu;
    sim_np = np;
    
    try
        % Save .mat files with the matrices
        % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
        'sim_Aobs','sim_Bobs','sim_Cobs','sim_Dobs','sim_TsObs','sim_Ts',...
        'sim_H','sim_K','sim_Hobs','sim_Kobs','sim_Cobs_fake','sim_Dobs_fake',...
        'sim_nx','sim_nu','sim_np','sim_ny','sim_nd','sim_ctrl',...
        'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
        
    catch
        error('Model cannot be updated. Maybe you are in the wrong folder.');
    end
    
    disp('Model succesfully updated.')
    disp(' ')
    disp('The model will now be opened...')
    disp(' ')
    cd(options.folder)
    pwaSys_model
    
else
    
    disp('Canceled.')
    disp(' ')
    
end






