function simulinkClosedLoopSim(object,options)
% simulinkClosedLoopSim       Generates the Simulink model of the 
%                             closed-loop dynamical system
%
% This is a private method called by method ltiSys/generateSimulinkModel.
%
% See also: ltiSys/generateSimulinkModel.


simulateVHDL = options.simulateVHDL;
circuit_parameters = options.circuit_parameters;
% Default value
choice = 'Overwrite';

% If the model already exists, ask what to do...
if exist([options.folder,'ltiSys_model.slx'],'file')
    
    % Construct a questdlg with two options
    choice = questdlg(['A model already exists in ',options.folder,'. Do you want to overwrite it?'], ...
        ' ', ...
        'Overwrite','Cancel','Cancel');
    
end

% Extract matrices
[A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy] = object.getMatrices();

% Extract controller
ctrl = object.controller;

% Number of variables
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;
nxref = numel(object.getController.getTrackingVariable);
trackVar = object.getController.getTrackingVariable;

% Put together matrices
sim_A = A; %#ok<*NASGU>
sim_B = [B Ex Fx Gx];
sim_C = [eye(object.nx);C];
sim_D = [zeros(object.nx,object.nu+object.np+object.nd+1);[D Ey Fy Gy]];

sim_xnames = object.getStateNames();
sim_unames = object.getInputNames();
sim_pnames = object.getParameterNames();
sim_ynames = object.getOutputNames();
sim_dnames = object.getUnmeasurableInputNames();

if strcmpi(choice,'Overwrite')    
    
    % Create folder
    mkdir(options.folder)
    
    cd(options.folder)
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_np = np;
    sim_nu = nu;
    sim_nxref = nxref;
    sim_ctrl = ctrl;
   
    % Copy S-function
    copyfile([options.simulinkFolder,'ct/s_MPCctrl.m'],[options.folder,'s_MPCctrl.m']);
    
    % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
        'sim_nx','sim_ny','sim_np','sim_nd','sim_nxref','sim_ctrl',...
        'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
    
    % Copy plot function
    copyfile([options.simulinkFolder,'ct/','plotScope.m'],[options.folder,'plotScope.m']);    
    
    % Generate a new Simulink model
    model = 'ltiSys_model';
    bdclose(model);
    new_system(model);

    %% LTI system generation
    add_block('simulink/Ports & Subsystems/Subsystem', [model, '/LTI System']);
    set_param([model, '/LTI System'], 'ContentPreviewEnabled', 'off');
    delete_line([model, '/LTI System'], 'In1/1', 'Out1/1');
    if object.isContinuousTime
        add_block('simulink/Continuous/State-Space', [model, '/LTI System/State-Space System']);
    else
        add_block('simulink/Discrete/Discrete State-Space', [model, '/LTI System/State-Space System']);
    end
    set_param([model, '/LTI System/State-Space System'], 'A', mat2str(sim_A), 'B', mat2str(sim_B), 'C', mat2str(sim_C), 'D', mat2str(sim_D), 'X0', 'sim_x0');
        
    if nxref ~= 0  
        add_block('simulink/Sources/Constant', [model, '/Constant reference']);
        set_param([model, '/Constant reference'], 'Value', 'sim_xref');
        set_param([model, '/Constant reference'], 'SampleTime', '-1');
        
        add_block('simulink/Sources/Signal Builder', [model, '/Reference builder']);
                
        add_block('simulink/Signal Routing/Manual Switch', [model, '/Reference switch']);
        add_line(model, 'Constant reference/1', 'Reference switch/1', 'autorouting', 'smart');
        add_line(model, 'Reference builder/1', 'Reference switch/2', 'autorouting', 'smart');
    end
    if np ~= 0  
        add_block('simulink/Sources/Constant', [model, '/Constant parameter']);
        set_param([model, '/Constant parameter'], 'Value', 'sim_p');
        set_param([model, '/Constant parameter'], 'SampleTime', '-1');
        
        add_block('simulink/Sources/Signal Builder', [model, '/Parameter builder']);
                
        add_block('simulink/Signal Routing/Manual Switch', [model, '/Parameter switch']);
        add_line(model, 'Constant parameter/1', 'Parameter switch/1', 'autorouting', 'smart');
        add_line(model, 'Parameter builder/1', 'Parameter switch/2', 'autorouting', 'smart');
    end
    if nd ~= 0  
        add_block('simulink/Sources/Constant', [model, '/Constant unmeasurable inputs']);
        set_param([model, '/Constant unmeasurable inputs'], 'Value', 'sim_d');
        set_param([model, '/Constant unmeasurable inputs'], 'SampleTime', '-1');
        
        add_block('simulink/Sources/Signal Builder', [model, '/Unmeasurable inputs builder']);
                
        add_block('simulink/Signal Routing/Manual Switch', [model, '/Unmeasurable inputs switch']);
        add_line(model, 'Constant unmeasurable inputs/1', 'Unmeasurable inputs switch/1', 'autorouting', 'smart');
        add_line(model, 'Unmeasurable inputs builder/1', 'Unmeasurable inputs switch/2', 'autorouting', 'smart');
    end
    
    % Define inputs and outputs of the system model
    if nu~=0
        set_param([model, '/LTI System/In1'], 'Name', 'u');
        add_block('simulink/Sources/Constant', [model, '/LTI System/Const']);
        add_block('simulink/Signal Routing/Mux', [model, '/LTI System/Mux']);
        if np~=0
            add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/p']);
            if nd~=0
                add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/d']);
                set_param([model, '/LTI System/Mux'], 'Inputs', '4');
                add_line([model, '/LTI System'], 'u/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'd/1', 'Mux/3', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Const/1', 'Mux/4', 'autorouting', 'smart');
                add_line(model, 'Parameter switch/1', 'LTI System/2', 'autorouting', 'smart');
                add_line(model, 'Unmeasurable inputs switch/1', 'LTI System/3', 'autorouting', 'smart');
            else
                set_param([model, '/LTI System/Mux'], 'Inputs', '3');
                add_line([model, '/LTI System'], 'u/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Const/1', 'Mux/3', 'autorouting', 'smart');
                add_line(model, 'Parameter switch/1', 'LTI System/2', 'autorouting', 'smart');
            end
        elseif nd~=0
            add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/d']);
            set_param([model, '/LTI System/Mux'], 'Inputs', '3');
            add_line([model, '/LTI System'], 'u/1', 'Mux/1', 'autorouting', 'smart');
            add_line([model, '/LTI System'], 'd/1', 'Mux/2', 'autorouting', 'smart');
            add_line([model, '/LTI System'], 'Const/1', 'Mux/3', 'autorouting', 'smart');
            add_line(model, 'Unmeasurable inputs switch/1', 'LTI System/2', 'autorouting', 'smart');
        else
            add_line([model, '/LTI System'], 'u/1', 'Mux/1', 'autorouting', 'smart');
            add_line([model, '/LTI System'], 'Const/1', 'Mux/2', 'autorouting', 'smart');
        end
        add_line([model, '/LTI System'], 'Mux/1', 'State-Space System/1', 'autorouting', 'smart');
    else
        delete_block([model, '/LTI System/In1']);
        add_block('simulink/Sources/Constant', [model, '/LTI System/Const']);
        if np~=0
            add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/p']);
            add_block('simulink/Signal Routing/Mux', [model, '/LTI System/Mux']);
            if nd~=0
                add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/d']);
                set_param([model, '/LTI System/Mux'], 'Inputs', '3');
                add_line([model, '/LTI System'], 'p/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'd/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Const/1', 'Mux/3', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Mux/1', 'State-Space System/1', 'autorouting', 'smart');               
                add_line(model, 'Parameter switch/1', 'LTI System/1', 'autorouting', 'smart');
                add_line(model, 'Unmeasurable inputs switch/1', 'LTI System/2', 'autorouting', 'smart');
            else
                add_line([model, '/LTI System'], 'p/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Const/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/LTI System'], 'Mux/1', 'State-Space System/1', 'autorouting', 'smart');
                add_line(model, 'Parameter switch/1', 'LTI System/1', 'autorouting', 'smart');
            end
        elseif nd~=0
            add_block('simulink/Ports & Subsystems/In1', [model, '/LTI System/d']);
            add_block('simulink/Signal Routing/Mux', [model, '/LTI System/Mux']);
            add_line([model, '/LTI System'], 'd/1', 'Mux/1', 'autorouting', 'smart');
            add_line([model, '/LTI System'], 'Const/1', 'Mux/2', 'autorouting', 'smart');
            add_line([model, '/LTI System'], 'Mux/1', 'State-Space System/1', 'autorouting', 'smart');
            add_line(model, 'Unmeasurable inputs switch/1', 'LTI System/1', 'autorouting', 'smart');
        else
            add_line([model, '/LTI System'], 'Const/1', 'State-Space System/1', 'autorouting', 'smart');
        end
    end
    
    if np ~= 0
        add_block('simulink/Sinks/Scope', [model, '/Monitor p']);
        add_line(model, 'Parameter switch/1', 'Monitor p/1', 'autorouting', 'smart');
        set_param([model, '/Monitor p'], 'BackgroundColor', 'Cyan');
        set_param([model, '/Monitor p'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_P', 'DataLoggingSaveFormat', 'Structure With Time');
    end
    
    if nd ~= 0
        add_block('simulink/Sinks/Scope', [model, '/Monitor d']);
        add_line(model, 'Unmeasurable inputs switch/1', 'Monitor d/1', 'autorouting', 'smart');
        set_param([model, '/Monitor d'], 'BackgroundColor', 'Cyan');
        set_param([model, '/Monitor d'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_D', 'DataLoggingSaveFormat', 'Structure With Time');
    end

    set_param([model, '/LTI System/Out1'], 'Name', 'x');
    add_block('simulink/Ports & Subsystems/Out1', [model, '/LTI System/y']);
    add_block('simulink/Signal Routing/Demux', [model, '/LTI System/Demux']);
    add_line([model, '/LTI System'], 'State-Space System/1', 'Demux/1', 'autorouting', 'smart');
    add_line([model, '/LTI System'], 'Demux/1', 'x/1', 'autorouting', 'smart');
    add_line([model, '/LTI System'], 'Demux/2', 'y/1', 'autorouting', 'smart');

    % Arrange the model and add scopes for u, x, y
    Simulink.BlockDiagram.arrangeSystem([model, '/LTI System']);
    add_block('simulink/Sinks/Scope', [model, '/Monitor x']);
    set_param([model, '/Monitor x'], 'BackgroundColor', 'Cyan');
    set_param([model, '/Monitor x'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_X', 'DataLoggingSaveFormat', 'Structure With Time');
    add_block('simulink/Sinks/Scope', [model, '/Monitor y']);
    set_param([model, '/Monitor y'], 'BackgroundColor', 'Cyan');
    set_param([model, '/Monitor y'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Y', 'DataLoggingSaveFormat', 'Structure With Time');
    add_line(model, 'LTI System/1', 'Monitor x/1', 'autorouting', 'smart');
    add_line(model, 'LTI System/2', 'Monitor y/1', 'autorouting', 'smart');
    
    
    %% Controller generation
    if object.hasController()

        % reset the lastwarn message and id
        lastwarn('', '');
        % try generating a System Generator token
        if options.simulateVHDL           
            add_block('sysgenController/ System Generator', [model, '/ System Generator']);
        end
        % check if a warning was raised
        [~, warnId] = lastwarn();
        if ~isempty(warnId)
            options.simulateVHDL = 0;
            warndlg('Xilinx System Generator is not properly installed. A Simulink model will be generated.');
        end

        % Generate a Simulink subsystem for the controller
        if ~options.simulateVHDL
            add_block('simulink/Ports & Subsystems/Subsystem', [model, '/Controller']);
            set_param([model, '/Controller'], 'ContentPreviewEnabled', 'off');
            delete_line([model, '/Controller'], 'In1/1', 'Out1/1');
            set_param([model, '/Controller/In1'], 'Name', 'x');
            set_param([model, '/Controller/Out1'], 'Name', 'u');
            add_block('simulink/User-Defined Functions/S-Function', [model, '/Controller/S-Function']);
            set_param([model, '/Controller/S-Function'], 'functionName', 's_MPCctrl', 'Parameters', 'sim_ctrl');
            if sim_ctrl.isTracking               
                add_block('simulink/Signal Routing/Mux', [model, '/Controller/Mux']);
                add_line([model, '/Controller'], 'x/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/Controller'], 'Mux/1', 'S-Function/1', 'autorouting', 'smart');
                if (np ~= 0) && (nd ~= 0)
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/p']);
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/d']);
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/xref']);
                    add_line(model, 'Parameter switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line(model, 'Unmeasurable inputs switch/1', 'Controller/3', 'autorouting', 'smart');
                    add_line(model, 'Reference switch/1', 'Controller/4', 'autorouting', 'smart');
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_np sim_nd sim_nxref]');
                    add_line([model, '/Controller'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'd/1', 'Mux/3', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'xref/1', 'Mux/4', 'autorouting', 'smart'); 
                elseif (np ~= 0)
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/p']);
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/xref']);
                    add_line(model, 'Parameter switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line(model, 'Reference switch/1', 'Controller/3', 'autorouting', 'smart');
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_np sim_nxref]');
                    add_line([model, '/Controller'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'xref/1', 'Mux/3', 'autorouting', 'smart');
                elseif (nd ~= 0)
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/d']);
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/xref']);
                    add_line(model, 'Unmeasurable inputs switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line(model, 'Reference switch/1', 'Controller/3', 'autorouting', 'smart');
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_nd sim_nxref]');
                    add_line([model, '/Controller'], 'd/1', 'Mux/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'xref/1', 'Mux/3', 'autorouting', 'smart');                    
                else
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/xref']);
                    add_line(model, 'Reference switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'xref/1', 'Mux/2', 'autorouting', 'smart');
                end
            else
                if (np ~= 0) && (nd ~= 0)
                    add_block('simulink/Signal Routing/Mux', [model, '/Controller/Mux']);
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_np sim_nd]');
                    add_line([model, '/Controller'], 'Mux/1', 'S-Function/1', 'autorouting', 'smart');
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/p']);
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/d']);
                    add_line(model, 'Parameter switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line(model, 'Unmeasurable inputs switch/1', 'Controller/3', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'x/1', 'Mux/1', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'd/1', 'Mux/3', 'autorouting', 'smart');
                elseif (np ~= 0)
                    add_block('simulink/Signal Routing/Mux', [model, '/Controller/Mux']);
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_np]');
                    add_line([model, '/Controller'], 'Mux/1', 'S-Function/1', 'autorouting', 'smart');
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/p']);
                    add_line(model, 'Parameter switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'x/1', 'Mux/1', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                elseif (nd ~= 0)
                    add_block('simulink/Signal Routing/Mux', [model, '/Controller/Mux']);
                    set_param([model, '/Controller/Mux'], 'Inputs', '[sim_nx sim_nd]');
                    add_line([model, '/Controller'], 'Mux/1', 'S-Function/1', 'autorouting', 'smart');
                    add_block('simulink/Ports & Subsystems/In1', [model, '/Controller/d']);
                    add_line(model, 'Unmeasurable inputs switch/1', 'Controller/2', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'x/1', 'Mux/1', 'autorouting', 'smart');
                    add_line([model, '/Controller'], 'd/1', 'Mux/2', 'autorouting', 'smart');                  
                else
                    add_line([model, '/Controller'], 'x/1', 'S-Function/1', 'autorouting', 'smart');
                end
            end
            add_line([model, '/Controller'], 'S-Function/1', 'u/1', 'autorouting', 'smart');
            add_line(model, 'LTI System/1', 'Controller/1', 'autorouting', 'smart');
            add_line(model, 'Controller/1', 'LTI System/1', 'autorouting', 'smart');
            add_block('simulink/Sinks/Scope', [model, '/Monitor u']);
            add_line(model, 'Controller/1', 'Monitor u/1', 'autorouting', 'smart');
            set_param([model, '/Monitor u'], 'BackgroundColor', 'Cyan');
            set_param([model, '/Monitor u'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_U', 'DataLoggingSaveFormat', 'Structure With Time');
            
            if sim_ctrl.isTracking
                add_block('simulink/Sinks/Scope', [model, '/Monitor xref']);
                add_line(model, 'Reference switch/1', 'Monitor xref/1', 'autorouting', 'smart');
                set_param([model, '/Monitor xref'], 'BackgroundColor', 'Cyan');
                set_param([model, '/Monitor xref'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Xref', 'DataLoggingSaveFormat', 'Structure With Time');
            end
            
            Simulink.BlockDiagram.arrangeSystem([model, '/Controller']);
                
        % Generate a Xilinx System Generator black box for the controller
        else
            
            % Default circuit parameters
            circuit_parameters = MPCctrlVHDLset(ctrl,options.circuit_parameters);           
            
            synthInfo = ctrl.generateVHDL(circuit_parameters);
            range = synthInfo.range;
            options.range = range;
            
            % Generate a configuration file for the black box
            writeSysgenConfigFile(object, options);
            
            % Frequency of the FPGA
            frequency = circuit_parameters.frequency;
            
            % Clock of the FPGA
            clk_period = 1/frequency;
            
            % System sampling time
            sysSamplingTime = ctrl.getInformation.sys.getSamplingTime;
            
            % Number of bits for the input (x)
            sim_n_bit_in = circuit_parameters.inputResolution;
            
            % Number of bits for the output (y)
            sim_n_bit_out = circuit_parameters.outputResolution;
            
            if circuit_parameters.useADC == 1
                sim_bin_point_in = 0;
            else
                sim_bin_point_in = sim_n_bit_in - ceil(log2(max(abs([range.xmax, range.xmin]),[],'all'))) - 1;
            end

            if circuit_parameters.useDAC == 1
                sim_bin_point_out = 0;
            else
                sim_bin_point_out = sim_n_bit_out - ceil(log2(max(abs([range.umax, range.umin]),[],'all'))) - 1;
            end

            % Compute arrays to transform the inputs from their actual (model) range to
            % the circuit range: x_cir = sim_x_scale_gain.*x + sim_x_scale_bias
            sim_x_scale_gain = (circuit_parameters.inputRange.max(1)-circuit_parameters.inputRange.min(1))./...
                (range.xmax(:)-range.xmin(:));
            sim_x_scale_bias = (range.xmax(:)+range.xmin(:))/2-...
                ((circuit_parameters.inputRange.max(1)+circuit_parameters.inputRange.min(1))/2./sim_x_scale_gain);

            if np == 0
                sim_p_scale_gain = [];
                sim_p_scale_bias = [];
            else
                sim_p_scale_gain = (synthInfo.circuit_parameters.inputRange.max(1)-synthInfo.circuit_parameters.inputRange.min(1))./...
                    (range.pmax(:)-range.pmin(:));
                sim_p_scale_bias = (range.pmax+range.pmin)/2-...
                    ((synthInfo.circuit_parameters.inputRange.max(1)+synthInfo.circuit_parameters.inputRange.min(1))/2./sim_p_scale_gain)';
            end
 
            if nd == 0
                sim_d_scale_gain = [];
                sim_d_scale_bias = [];
            else
                sim_d_scale_gain = (synthInfo.circuit_parameters.inputRange.max(1)-synthInfo.circuit_parameters.inputRange.min(1))./...
                    (range.dmax(:)-range.dmin(:));
                sim_d_scale_bias = (range.dmax+range.dmin)/2-...
                    ((synthInfo.circuit_parameters.inputRange.max(1)+synthInfo.circuit_parameters.inputRange.min(1))/2./sim_d_scale_gain)';
            end

            if nxref == 0
                sim_xref_scale_gain = 0;
                sim_xref_scale_bias = 0;
            else
                sim_xref_scale_gain = (synthInfo.circuit_parameters.inputRange.max(1)-synthInfo.circuit_parameters.inputRange.min(1))./...
                    (range.xrefmax(:)-range.xrefmin(:));
                sim_xref_scale_bias = (range.xrefmax(:)+range.xrefmin(:))/2-...
                    ((synthInfo.circuit_parameters.inputRange.max(1)+synthInfo.circuit_parameters.inputRange.min(1))/2./sim_xref_scale_gain);
            end

            % Compute arrays to transform the outputs from their circuit range to the
            % actual (model) range: u = sim_u_scale_gain.*u + sim_u_scale_bias
            sim_u_scale_gain = (range.umax(:)-range.umin(:))./(circuit_parameters.outputRange.max(1)-circuit_parameters.outputRange.min(1));
            sim_u_scale_bias = (range.umax(:)+range.umin(:))/2-...
                ((circuit_parameters.outputRange.max(1)+circuit_parameters.outputRange.min(1))/2.*sim_u_scale_gain)';

            % Save workspace for Simulink simulations
            save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
                'sim_nx','sim_ny','sim_np','sim_nd','sim_nxref',...
                'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames','clk_period',...
                'sim_n_bit_in','sim_n_bit_out','sim_nu','sim_bin_point_in','sim_bin_point_out',...
                'sim_u_scale_gain','sim_u_scale_bias','sim_x_scale_gain','sim_x_scale_bias',...
                'sim_p_scale_gain','sim_p_scale_bias','sim_d_scale_gain','sim_d_scale_bias',...
                'sim_xref_scale_gain','sim_xref_scale_bias');
            
%             % Add System Generator token
%             add_block('sysgenController/ System Generator', [model, '/ System Generator']);
            
            % Add black box
            add_block('sysgenController/Controller', [model, '/Controller']);

            % Add monitor for output (u)
            add_block('simulink/Sinks/Scope', [model, '/Monitor u']);
            set_param([model, '/Monitor u'], 'BackgroundColor', 'Cyan');
            set_param([model, '/Monitor u'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_U', 'DataLoggingSaveFormat', 'Structure With Time');

            % Add reset input (step)
            add_block('sysgenController/bool_input', [model, '/reset']);
            add_block('simulink/Sources/Step', [model, '/rst']);
            set_param([model, '/rst'], 'Time', 'clk_period/2', 'Before', '0', 'After', '1');
            add_line(model, 'rst/1', 'reset/1', 'autorouting', 'smart');
            add_line(model, 'reset/1', 'Controller/1', 'autorouting', 'smart');
            
            % Add start input (pulse)
            if isa(ctrl, 'implicitMPCctrl')
                add_block('sysgenController/bool_input', [model, '/start']);
                add_block('simulink/Sources/Pulse Generator', [model, '/str']);
                set_param([model, '/str'], 'Period', num2str(sysSamplingTime), 'PulseWidth', [num2str(100/sysSamplingTime), '*clk_period'], 'PhaseDelay', 'clk_period/2');
                add_line(model, 'str/1', 'start/1', 'autorouting', 'smart');
                add_line(model, 'start/1', 'Controller/2', 'autorouting', 'smart');
            end
            
            % Scale input (x)
            if circuit_parameters.useADC == 1
                add_block('simulink/Ports & Subsystems/Subsystem', [model, '/stateScale']);
                set_param([model, '/stateScale'], 'ContentPreviewEnabled', 'off');
                delete_line([model, '/stateScale'], 'In1/1', 'Out1/1');
                set_param([model, '/stateScale/In1'], 'Name', 'x');
                set_param([model, '/stateScale/Out1'], 'Name', 'x_cir');
                add_block('simulink/Math Operations/Add', [model, '/stateScale/Add']);
                set_param([model, '/stateScale/Add'], 'Inputs', '+-');
                add_block('simulink/Math Operations/Product', [model, '/stateScale/Product']);
                add_block('simulink/Sources/Constant', [model, '/stateScale/bias']);
                set_param([model, '/stateScale/bias'], 'Value', 'sim_x_scale_bias');
                add_block('simulink/Sources/Constant', [model, '/stateScale/gain']);
                set_param([model, '/stateScale/gain'], 'Value', 'sim_x_scale_gain');
                add_line([model, '/stateScale'], 'x/1', 'Add/1', 'autorouting', 'smart');
                add_line([model, '/stateScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                add_line([model, '/stateScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
                add_line([model, '/stateScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                add_line([model, '/stateScale'], 'Product/1', 'x_cir/1', 'autorouting', 'smart');
                Simulink.BlockDiagram.arrangeSystem([model, '/stateScale']);
                add_line(model, 'LTI System/1', 'stateScale/1', 'autorouting', 'smart');
            end
            
            % Scale parameters (p)
            if np ~= 0  
                if circuit_parameters.useADC == 1
                    add_block('simulink/Ports & Subsystems/Subsystem', [model, '/parScale']);
                    set_param([model, '/parScale'], 'ContentPreviewEnabled', 'off');
                    delete_line([model, '/parScale'], 'In1/1', 'Out1/1');
                    set_param([model, '/parScale/In1'], 'Name', 'p');
                    set_param([model, '/parScale/Out1'], 'Name', 'p_cir');
                    add_block('simulink/Math Operations/Add', [model, '/parScale/Add']);
                    set_param([model, '/parScale/Add'], 'Inputs', '+-');
                    add_block('simulink/Math Operations/Product', [model, '/parScale/Product']);
                    add_block('simulink/Sources/Constant', [model, '/parScale/bias']);
                    set_param([model, '/parScale/bias'], 'Value', 'sim_p_scale_bias');
                    add_block('simulink/Sources/Constant', [model, '/parScale/gain']);
                    set_param([model, '/parScale/gain'], 'Value', 'sim_p_scale_gain');
                    add_line([model, '/parScale'], 'p/1', 'Add/1', 'autorouting', 'smart');
                    add_line([model, '/parScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                    add_line([model, '/parScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
                    add_line([model, '/parScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                    add_line([model, '/parScale'], 'Product/1', 'p_cir/1', 'autorouting', 'smart');
                    Simulink.BlockDiagram.arrangeSystem([model, '/parScale']);
                    add_line(model, 'Parameter switch/1', 'parScale/1', 'autorouting', 'smart');
                end
            end
            
            % Scale unmeasurable inputs or parameters (d)
            if nd ~= 0  
                if circuit_parameters.useADC == 1
                    add_block('simulink/Ports & Subsystems/Subsystem', [model, '/unmeasInputScale']);
                    set_param([model, '/unmeasInputScale'], 'ContentPreviewEnabled', 'off');
                    delete_line([model, '/unmeasInputScale'], 'In1/1', 'Out1/1');
                    set_param([model, '/unmeasInputScale/In1'], 'Name', 'd');
                    set_param([model, '/unmeasInputScale/Out1'], 'Name', 'd_cir');
                    add_block('simulink/Math Operations/Add', [model, '/unmeasInputScale/Add']);
                    set_param([model, '/unmeasInputScale/Add'], 'Inputs', '+-');
                    add_block('simulink/Math Operations/Product', [model, '/unmeasInputScale/Product']);
                    add_block('simulink/Sources/Constant', [model, '/unmeasInputScale/bias']);
                    set_param([model, '/unmeasInputScale/bias'], 'Value', 'sim_d_scale_bias');
                    add_block('simulink/Sources/Constant', [model, '/unmeasInputScale/gain']);
                    set_param([model, '/unmeasInputScale/gain'], 'Value', 'sim_d_scale_gain');
                    add_line([model, '/unmeasInputScale'], 'd/1', 'Add/1', 'autorouting', 'smart');
                    add_line([model, '/unmeasInputScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                    add_line([model, '/unmeasInputScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
                    add_line([model, '/unmeasInputScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                    add_line([model, '/unmeasInputScale'], 'Product/1', 'd_cir/1', 'autorouting', 'smart');
                    Simulink.BlockDiagram.arrangeSystem([model, '/unmeasInputScale']);
                    add_line(model, 'Unmeasurable inputs switch/1', 'unmeasInputScale/1', 'autorouting', 'smart');
                end
            end
            
            % Scale reference input (xref)
            if nxref ~= 0  
                if circuit_parameters.useADC == 1
                    add_block('simulink/Ports & Subsystems/Subsystem', [model, '/refScale']);
                    set_param([model, '/refScale'], 'ContentPreviewEnabled', 'off');
                    delete_line([model, '/refScale'], 'In1/1', 'Out1/1');
                    set_param([model, '/refScale/In1'], 'Name', 'xref');
                    set_param([model, '/refScale/Out1'], 'Name', 'xref_cir');
                    add_block('simulink/Math Operations/Add', [model, '/refScale/Add']);
                    set_param([model, '/refScale/Add'], 'Inputs', '+-');
                    add_block('simulink/Math Operations/Product', [model, '/refScale/Product']);
                    add_block('simulink/Sources/Constant', [model, '/refScale/bias']);
                    set_param([model, '/refScale/bias'], 'Value', 'sim_xref_scale_bias');
                    add_block('simulink/Sources/Constant', [model, '/refScale/gain']);
                    set_param([model, '/refScale/gain'], 'Value', 'sim_xref_scale_gain');
                    add_line([model, '/refScale'], 'xref/1', 'Add/1', 'autorouting', 'smart');
                    add_line([model, '/refScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                    add_line([model, '/refScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
                    add_line([model, '/refScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                    add_line([model, '/refScale'], 'Product/1', 'xref_cir/1', 'autorouting', 'smart');
                    Simulink.BlockDiagram.arrangeSystem([model, '/refScale']);
                    add_line(model, 'Reference switch/1', 'refScale/1', 'autorouting', 'smart');
                end
            end
                
            % Scale output (u)
            if circuit_parameters.useDAC == 1
                add_block('simulink/Ports & Subsystems/Subsystem', [model, '/controlScale']);
                set_param([model, '/controlScale'], 'ContentPreviewEnabled', 'off');
                delete_line([model, '/controlScale'], 'In1/1', 'Out1/1');
                set_param([model, '/controlScale/In1'], 'Name', 'u_cir');
                set_param([model, '/controlScale/Out1'], 'Name', 'u');
                add_block('simulink/Math Operations/Product', [model, '/controlScale/Product']);
                add_block('simulink/Math Operations/Add', [model, '/controlScale/Add']);
                add_block('simulink/Sources/Constant', [model, '/controlScale/gain']);
                set_param([model, '/controlScale/gain'], 'Value', 'sim_u_scale_gain');
                add_block('simulink/Sources/Constant', [model, '/controlScale/bias']);
                set_param([model, '/controlScale/bias'], 'Value', 'sim_u_scale_bias');
                add_line([model, '/controlScale'], 'u_cir/1', 'Product/1', 'autorouting', 'smart');
                add_line([model, '/controlScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                add_line([model, '/controlScale'], 'Product/1', 'Add/1', 'autorouting', 'smart');
                add_line([model, '/controlScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                add_line([model, '/controlScale'], 'Add/1', 'u/1', 'autorouting', 'smart');
                Simulink.BlockDiagram.arrangeSystem([model, '/controlScale']);
                add_line(model, 'controlScale/1', 'LTI System/1', 'autorouting', 'smart');
            end
               
            % Define inputs and outputs of the controller
            if nx == 1
                if strcmp(circuit_parameters.inputRepresentation,'signed')
                    add_block('sysgenController/signed_fixed_input', [model, '/x']);
                else
                    add_block('sysgenController/unsigned_fixed_input', [model, '/x']);
                end
                if circuit_parameters.useADC == 1
                    add_line(model, 'stateScale/1', 'x/1', 'autorouting', 'smart');
                else
                    add_line(model, 'LTI System/1', 'x/1', 'autorouting', 'smart');
                end
                if isa(ctrl, 'implicitMPCctrl')
                    add_line(model, 'x/1', 'Controller/3', 'autorouting', 'smart');
                else
                    add_line(model, 'x/1', 'Controller/2', 'autorouting', 'smart');
                end
            else
                add_block('simulink/Signal Routing/Demux', [model, '/Demux']);
                set_param([model, '/Demux'], 'Outputs', num2str(nx));
                if circuit_parameters.useADC == 1
                    add_line(model, 'stateScale/1', 'Demux/1', 'autorouting', 'smart');
                else
                    add_line(model, 'LTI System/1', 'Demux/1', 'autorouting', 'smart');
                end
                for i=1:nx
                    if strcmp(circuit_parameters.inputRepresentation,'signed')
                        add_block('sysgenController/signed_fixed_input', [model, '/x', num2str(i)]);
                    else
                        add_block('sysgenController/unsigned_fixed_input', [model, '/x', num2str(i)]);
                    end
                    add_line(model, ['Demux/', num2str(i)], ['x', num2str(i), '/1'], 'autorouting', 'smart');
                    if isa(ctrl, 'implicitMPCctrl')
                        add_line(model, ['x', num2str(i), '/1'], ['Controller/', num2str(i+2)], 'autorouting', 'smart');
                    else
                        add_line(model, ['x', num2str(i), '/1'], ['Controller/', num2str(i+1)], 'autorouting', 'smart');
                    end
                end
            end
            
            if np ~= 0
                if np == 1
                    if strcmp(circuit_parameters.inputRepresentation,'signed')
                        add_block('sysgenController/signed_fixed_input', [model, '/p']);
                    else
                        add_block('sysgenController/unsigned_fixed_input', [model, '/p']);
                    end
                    if circuit_parameters.useADC == 1
                        add_line(model, 'parScale/1', 'p/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Parameter switch/1', 'p/1', 'autorouting', 'smart');
                    end
                    if isa(ctrl, 'implicitMPCctrl')
                        add_line(model, 'p/1', ['Controller/', num2str(nx+3)], 'autorouting', 'smart');
                    else
                        add_line(model, 'p/1', ['Controller/', num2str(nx+2)], 'autorouting', 'smart');
                    end
                else
                    add_block('simulink/Signal Routing/Demux', [model, '/DemuxP']);
                    set_param([model, '/DemuxP'], 'Outputs', num2str(np));
                    if circuit_parameters.useADC == 1
                        add_line(model, 'parScale/1', 'DemuxP/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Parameter switch/1', 'DemuxP/1', 'autorouting', 'smart');
                    end
                    for i=1:np         
                        if strcmp(circuit_parameters.inputRepresentation,'signed')
                            add_block('sysgenController/signed_fixed_input', [model, '/p', num2str(i)]);
                        else
                            add_block('sysgenController/unsigned_fixed_input', [model, '/p', num2str(i)]);
                        end
                        add_line(model, ['DemuxP/', num2str(i)], ['p', num2str(i), '/1'], 'autorouting', 'smart');
                        if isa(ctrl, 'implicitMPCctrl')
                            add_line(model, ['p', num2str(i), '/1'], ['Controller/', num2str(nx+2+i)], 'autorouting', 'smart');
                        else
                            add_line(model, ['p', num2str(i), '/1'], ['Controller/', num2str(nx+1+i)], 'autorouting', 'smart');
                        end
                    end
                end
            end
            
            if nd ~= 0
                if nd == 1
                    if strcmp(circuit_parameters.inputRepresentation,'signed')
                        add_block('sysgenController/signed_fixed_input', [model, '/d']);
                    else
                        add_block('sysgenController/unsigned_fixed_input', [model, '/d']);
                    end
                    if circuit_parameters.useADC == 1
                        add_line(model, 'unmeasInputScale/1', 'd/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Unmeasurable inputs switch/1', 'd/1', 'autorouting', 'smart');
                    end
                    if isa(ctrl, 'implicitMPCctrl')
                        add_line(model, 'd/1', ['Controller/', num2str(nx+np+3)], 'autorouting', 'smart');
                    else
                        add_line(model, 'd/1', ['Controller/', num2str(nx+np+2)], 'autorouting', 'smart');
                    end
                else
                    add_block('simulink/Signal Routing/Demux', [model, '/DemuxD']);
                    set_param([model, '/DemuxD'], 'Outputs', num2str(nd));
                    if circuit_parameters.useADC == 1
                        add_line(model, 'unmeasInputScale/1', 'DemuxD/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Unmeasurable inputs switch/1', 'DemuxD/1', 'autorouting', 'smart');
                    end
                    for i=1:nd
                        if strcmp(circuit_parameters.inputRepresentation,'signed')
                            add_block('sysgenController/signed_fixed_input', [model, '/d', num2str(i)]);
                        else
                            add_block('sysgenController/unsigned_fixed_input', [model, '/d', num2str(i)]);
                        end
                        add_line(model, ['DemuxD/', num2str(i)], ['d', num2str(i), '/1'], 'autorouting', 'smart');
                        if isa(ctrl, 'implicitMPCctrl')
                            add_line(model, ['d', num2str(i), '/1'], ['Controller/', num2str(nx+np+2+i)], 'autorouting', 'smart');
                        else
                            add_line(model, ['d', num2str(i), '/1'], ['Controller/', num2str(nx+np+1+i)], 'autorouting', 'smart');
                        end
                    end
                end
            end            
            
            if nxref ~= 0
                if nxref == 1
                    if strcmp(circuit_parameters.inputRepresentation,'signed')
                        add_block('sysgenController/signed_fixed_input', [model, '/xref']);
                    else
                        add_block('sysgenController/unsigned_fixed_input', [model, '/xref']);
                    end
                    if circuit_parameters.useADC == 1
                        add_line(model, 'refScale/1', 'xref/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Reference switch/1', 'xref/1', 'autorouting', 'smart');
                    end
                    if isa(ctrl, 'implicitMPCctrl')
                        add_line(model, 'xref/1', ['Controller/', num2str(nx+np+nd+3)], 'autorouting', 'smart');
                    else
                        add_line(model, 'xref/1', ['Controller/', num2str(nx+np+nd+2)], 'autorouting', 'smart');
                    end
                else
                    add_block('simulink/Signal Routing/Demux', [model, '/Demux1']);
                    set_param([model, '/Demux1'], 'Outputs', num2str(nxref));
                    if circuit_parameters.useADC == 1
                        add_line(model, 'refScale/1', 'Demux1/1', 'autorouting', 'smart');
                    else
                        add_line(model, 'Reference switch/1', 'Demux1/1', 'autorouting', 'smart');
                    end
                    for i=1:nxref 
                        if strcmp(circuit_parameters.inputRepresentation,'signed')        
                            add_block('sysgenController/signed_fixed_input', [model, '/xref', num2str(i)]);
                        else
                            add_block('sysgenController/unsigned_fixed_input', [model, '/xref', num2str(i)]);
                        end
                        add_line(model, ['Demux1/', num2str(i)], ['xref', num2str(i), '/1'], 'autorouting', 'smart');
                        if isa(ctrl, 'implicitMPCctrl')
                            add_line(model, ['xref', num2str(i), '/1'], ['Controller/', num2str(nx+np+nd+2+i)], 'autorouting', 'smart');
                        else
                            add_line(model, ['xref', num2str(i), '/1'], ['Controller/', num2str(nx+np+nd+1+i)], 'autorouting', 'smart');
                        end
                    end
                end
            end
            
            if nu == 1
                add_block('xbsBasic_r4/Gateway Out', [model, '/u']);
                add_line(model, 'Controller/1', 'u/1', 'autorouting', 'smart');
                if circuit_parameters.useDAC == 1
                    add_line(model, 'u/1', 'controlScale/1', 'autorouting', 'smart');
                    add_line(model, 'controlScale/1', 'Monitor u/1', 'autorouting', 'smart');
                else
                    add_line(model, 'u/1', 'LTI System/1', 'autorouting', 'smart');
                    add_line(model, 'u/1', 'Monitor u/1', 'autorouting', 'smart');
                end
            else
                add_block('simulink/Signal Routing/Mux', [model, '/Mux']);
                set_param([model, '/Mux'], 'Inputs', num2str(nu));
                if circuit_parameters.useDAC == 1
                    add_line(model, 'Mux/1', 'controlScale/1', 'autorouting', 'smart');
                    add_line(model, 'controlScale/1', 'Monitor u/1', 'autorouting', 'smart');
                else
                    add_line(model, 'Mux/1', 'LTI System/1', 'autorouting', 'smart');
                    add_line(model, 'Mux/1', 'Monitor u/1', 'autorouting', 'smart');
                end
                for i=1:nu
                    add_block('xbsBasic_r4/Gateway Out', [model, '/u', num2str(i)]);
                    add_line(model, ['Controller/', num2str(i)], ['u', num2str(i), '/1'], 'autorouting', 'smart');
                    add_line(model, ['u', num2str(i), '/1'], ['Mux/', num2str(i)], 'autorouting', 'smart');
                end
            end
            
            % Add outputs 'sample' and 'done'
            add_block('xbsBasic_r4/Gateway Out', [model, '/sample']);
            add_block('xbsBasic_r4/Gateway Out', [model, '/done']);
            add_block('simulink/Sinks/Scope', [model, '/Monitor sample']);
            add_block('simulink/Sinks/Scope', [model, '/Monitor done']);
            add_line(model, ['Controller/', num2str(nu+1)], 'sample/1', 'autorouting', 'smart');
            add_line(model, ['Controller/', num2str(nu+2)], 'done/1', 'autorouting', 'smart');
            add_line(model, 'sample/1', 'Monitor sample/1', 'autorouting', 'smart');
            add_line(model, 'done/1', 'Monitor done/1', 'autorouting', 'smart');
        end
    end
    
    % Arrange the model and zoom on it
    Simulink.BlockDiagram.arrangeSystem(model);
    set_param(model, 'ZoomFactor','FitSystem')

    % Hide subsystems inside
    set_param(getActiveConfigSet(model), 'ReturnWorkspaceOutputs', 'Off');
    set_param(getActiveConfigSet(model), 'StopTime', 'sim_T');
    
    % Define callback functions for the model
    set_param(model, 'PostLoadFcn', "load simulink.mat"+newline+"setParameters");
    set_param(model, 'InitFcn', "load simulink.mat"+newline+"setParameters");
    set_param(model, 'StartFcn', "load simulink.mat"+newline+"setParameters"+newline+"clear res_X res_U res_Y res_P res_D res_Xref"+newline+"clear res_Xest res_Dest res_Yest");
    set_param(model, 'StopFcn', "plotScope");

    % Save Simulink model
    save_system(model, [options.folder,'ltiSys_model.slx']);

    
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
    edit('setParameters.m')
    open_system(model);
    
    if ~options.simulateVHDL
        helpdlg('Before running the model, set the parameters in file "setParameters.m"')
    else
        helpdlg("Before running the model:" + newline + "- select the FPGA device in System Generator" + newline + "- set the parameters in file ""setParameters.m""");
    end
    
else
    
    disp('Canceled.')
    disp(' ')
    
end

end

