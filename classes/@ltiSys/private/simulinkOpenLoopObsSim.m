function simulinkOpenLoopObsSim(object,options)
% simulinkOpenLoopObsSim       Generates the Simulink model of the 
%                              open loop dynamical system with observer
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

% Extract observer
obs = object.getObserver();

% Extract observer matrices, gain and sampling time
[Aobs, Bobs, Cobs, Dobs, Gxobs, Gyobs] = obs.getMatrices();

% Number of variables
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;

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

sim_TsObs = obs.getSamplingTime();

% Put together observer matrices
sim_Aobs = Aobs;
sim_Bobs = [Bobs Gxobs];
sim_Cobs = Cobs;
sim_Dobs = [Dobs Gyobs];

if strcmpi(choice,'Overwrite')    
    
    preFolder = '';
    % Create folder
    mkdir(options.folder)
    
    cd(options.folder)
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_np = np;
    sim_nu = nu;
    
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
        end
        add_line([model, '/LTI System'], 'Mux/1', 'State-Space System/1', 'autorouting', 'smart');
        add_block('simulink/Sources/Constant', [model, '/Constant input']);
        set_param([model, '/Constant input'], 'Value', 'sim_u');
        
        add_block('simulink/Sources/Signal Builder', [model, '/Input builder']);
                
        add_block('simulink/Signal Routing/Manual Switch', [model, '/Input switch']);
        add_line(model, 'Constant input/1', 'Input switch/1', 'autorouting', 'smart');
        add_line(model, 'Input builder/1', 'Input switch/2', 'autorouting', 'smart');
        
        add_line(model, 'Input switch/1', 'LTI System/1', 'autorouting', 'smart');
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
    
    add_block('simulink/Sinks/Scope', [model, '/Monitor xest']);
    set_param([model, '/Monitor xest'], 'BackgroundColor', 'Cyan');
    set_param([model, '/Monitor xest'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Xest', 'DataLoggingSaveFormat', 'Structure With Time');
    
    %% Observer generation
    if ~simulateVHDL
        % Save .mat files with the matrices
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
            'sim_Aobs','sim_Bobs','sim_Cobs','sim_Dobs','sim_TsObs',...
            'sim_nx','sim_nu','sim_ny','sim_nd',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
        
        add_block('simulink/Ports & Subsystems/Subsystem', [model, '/Observer']);
        set_param([model, '/Observer'], 'ContentPreviewEnabled', 'off');
        delete_line([model, '/Observer'], 'In1/1', 'Out1/1');
        
        add_block('simulink/Discrete/Discrete State-Space', [model, '/Observer/Discrete State-Space']);
        set_param([model, '/Observer/Discrete State-Space'], 'A', 'sim_Aobs', 'B', 'sim_Bobs', 'C', 'sim_Cobs', 'D', 'sim_Dobs', 'X0', '[sim_x0obs; sim_d0obs]', 'SampleTime', 'sim_TsObs');
        
        % Define inputs and outputs of the observer
        
        set_param([model, '/Observer/In1'], 'Name', 'y');
        add_line(model, 'LTI System/2', 'Observer/1', 'autorouting', 'smart');
        add_block('simulink/Sources/Constant', [model, '/Observer/Const']);
        add_block('simulink/Signal Routing/Mux', [model, '/Observer/Mux']);
        if nu~=0
            add_block('simulink/Ports & Subsystems/In1', [model, '/Observer/u']);
            add_line(model, 'Constant input/1', 'Observer/2', 'autorouting', 'smart');
            if np~=0
                add_block('simulink/Ports & Subsystems/In1', [model, '/Observer/p']);
                set_param([model, '/Observer/Mux'], 'Inputs', '4');
                add_line([model, '/Observer'], 'u/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'p/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'y/1', 'Mux/3', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'Const/1', 'Mux/4', 'autorouting', 'smart');
                add_line(model, 'Parameter switch/1', 'Observer/3', 'autorouting', 'smart');
            else
                set_param([model, '/Observer/Mux'], 'Inputs', '3');
                add_line([model, '/Observer'], 'u/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'y/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'Const/1', 'Mux/3', 'autorouting', 'smart');
            end
        else
            if np~=0
                add_block('simulink/Ports & Subsystems/In1', [model, '/Observer/p']);
                set_param([model, '/Observer/Mux'], 'Inputs', '3');
                add_line([model, '/Observer'], 'p/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'y/1', 'Mux/2', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'Const/1', 'Mux/3', 'autorouting', 'smart');
                add_line(model, 'Parameter switch/1', 'Observer/2', 'autorouting', 'smart');
            else
                add_line([model, '/Observer'], 'y/1', 'Mux/1', 'autorouting', 'smart');
                add_line([model, '/Observer'], 'Const/1', 'Mux/2', 'autorouting', 'smart');
            end
        end
        add_line([model, '/Observer'], 'Mux/1', 'Discrete State-Space/1', 'autorouting', 'smart');
        set_param([model, '/Observer/Out1'], 'Name', 'xest');
        add_block('simulink/Signal Routing/Demux', [model, '/Observer/Demux']);
        if nd~=0
            set_param([model, '/Observer/Demux'], 'Outputs', '[sim_ny sim_nx sim_nd]');
        else
            set_param([model, '/Observer/Demux'], 'Outputs', '[sim_ny sim_nx]');
        end
        
        % Kalman predictor
        if isa(obs, 'kalmanPredictor')
            add_line([model, '/Observer'], 'Discrete State-Space/1', 'Demux/1', 'autorouting', 'smart');
            
        % Kalman filter
        else
            set_param([model, '/Observer/Discrete State-Space'], 'C', 'eye(size(sim_Aobs))', 'D', 'zeros(size(sim_Aobs,1),size(sim_Bobs,2))');
            add_block('simulink/Ports & Subsystems/Subsystem', [model, '/Observer/Compute outputs']);
            set_param([model, '/Observer/Compute outputs'], 'ContentPreviewEnabled', 'off');
            delete_line([model, '/Observer/Compute outputs'], 'In1/1', 'Out1/1');
            set_param([model, '/Observer/Compute outputs/In1'], 'Name', 'x');
            add_block('simulink/Ports & Subsystems/In1', [model, '/Observer/Compute outputs/in']);
            set_param([model, '/Observer/Compute outputs/Out1'], 'Name', 'out');
            add_block('simulink/Signal Routing/Demux', [model, '/Observer/Compute outputs/Demux']);
            add_block('simulink/Signal Routing/Mux', [model, '/Observer/Compute outputs/Mux']);
            add_line([model, '/Observer/Compute outputs'], 'in/1', 'Demux/1', 'autorouting', 'smart');
            if np~=0
                set_param([model, '/Observer/Compute outputs/Demux'], 'Outputs', '[sim_nu sim_np sim_ny 1]');
                set_param([model, '/Observer/Compute outputs/Mux'], 'Inputs', '3');
                add_line([model, '/Observer/Compute outputs'], 'Demux/4', 'Mux/3', 'autorouting', 'smart');
            else
                set_param([model, '/Observer/Compute outputs/Demux'], 'Outputs', '[sim_nu sim_ny 1]');
                set_param([model, '/Observer/Compute outputs/Mux'], 'Inputs', '2');
            end
            add_block('simulink/Sinks/Terminator', [model, '/Observer/Compute outputs/Terminator']);
            add_line([model, '/Observer/Compute outputs'], 'Demux/1', 'Terminator/1', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Demux/2', 'Mux/1', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Demux/3', 'Mux/2', 'autorouting', 'smart');
            add_block('simulink/Quick Insert/Math Operations/Matrix Multiply', [model, '/Observer/Compute outputs/Matrix multiply 1']);
            set_param([model, '/Observer/Compute outputs/Matrix multiply 1'], 'SampleTime', 'sim_TsObs');
            add_block('simulink/Quick Insert/Math Operations/Matrix Multiply', [model, '/Observer/Compute outputs/Matrix multiply 2']);
            set_param([model, '/Observer/Compute outputs/Matrix multiply 2'], 'SampleTime', 'sim_TsObs');
            add_block('simulink/Sources/Constant', [model, '/Observer/Compute outputs/Constant 1']);
            set_param([model, '/Observer/Compute outputs/Constant 1'], 'Value', 'sim_Cobs');
            add_block('simulink/Sources/Constant', [model, '/Observer/Compute outputs/Constant 2']);
            set_param([model, '/Observer/Compute outputs/Constant 2'], 'Value', 'sim_Dobs(:,sim_nu+1:end)');
            add_line([model, '/Observer/Compute outputs'], 'Constant 1/1', 'Matrix multiply 1/1', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'x/1', 'Matrix multiply 1/2', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Constant 2/1', 'Matrix multiply 2/1', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Mux/1', 'Matrix multiply 2/2', 'autorouting', 'smart');
            add_block('simulink/Math Operations/Sum', [model, '/Observer/Compute outputs/Sum']);
            add_line([model, '/Observer/Compute outputs'], 'Matrix multiply 1/1', 'Sum/1', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Matrix multiply 2/1', 'Sum/2', 'autorouting', 'smart');
            add_line([model, '/Observer/Compute outputs'], 'Sum/1', 'out/1', 'autorouting', 'smart');
            
            Simulink.BlockDiagram.arrangeSystem([model, '/Observer/Compute outputs']);
            
            add_line([model, '/Observer'], 'Discrete State-Space/1', 'Compute outputs/1', 'autorouting', 'smart');
            add_line([model, '/Observer'], 'Mux/1', 'Compute outputs/2', 'autorouting', 'smart');
            add_line([model, '/Observer'], 'Compute outputs/1', 'Demux/1', 'autorouting', 'smart');
        end
        
        add_block('simulink/Sinks/Scope', [model, '/Monitor yest']);
        set_param([model, '/Monitor yest'], 'BackgroundColor', 'Cyan');
        set_param([model, '/Monitor yest'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Yest', 'DataLoggingSaveFormat', 'Structure With Time');
        add_line([model, '/Observer'], 'Demux/2', 'xest/1', 'autorouting', 'smart');
        if nd~=0
            set_param([model, '/Observer/Demux'], 'Outputs', '3');
            add_block('simulink/Ports & Subsystems/Out1', [model, '/Observer/dest']);
            add_block('simulink/Ports & Subsystems/Out1', [model, '/Observer/yest']);
            add_line([model, '/Observer'], 'Demux/3', 'dest/1', 'autorouting', 'smart');
            add_block('simulink/Sinks/Scope', [model, '/Monitor dest']);
            set_param([model, '/Monitor dest'], 'BackgroundColor', 'Cyan');
            set_param([model, '/Monitor dest'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Dest', 'DataLoggingSaveFormat', 'Structure With Time');
            add_line(model, 'Observer/2', 'Monitor dest/1', 'autorouting', 'smart');
            add_line(model, 'Observer/3', 'Monitor yest/1', 'autorouting', 'smart');
        else
            add_block('simulink/Ports & Subsystems/Out1', [model, '/Observer/yest']);
            add_line(model, 'Observer/2', 'Monitor yest/1', 'autorouting', 'smart');            
        end
        add_line([model, '/Observer'], 'Demux/1', 'yest/1', 'autorouting', 'smart');
        
        Simulink.BlockDiagram.arrangeSystem([model, '/Observer']);
        add_line(model, 'Observer/1', 'Monitor xest/1', 'autorouting', 'smart');
    
    % Generate a Xilinx System Generator black box for the observer
    else
        
        % Default circuit parameters
        range = circuit_parameters.range;
        circuit_parameters = MPCctrlVHDLset(obs,options.circuit_parameters);
        circuit_parameters.range = range;
        synthInfo = obs.generateVHDL(circuit_parameters);

        options.range = range; 
        
        % Generate a configuration file for the black box
        writeObserverSysgenConfigFile(object, options);
        
        % Frequency of the FPGA
        frequency = circuit_parameters.frequency;

        % Clock of the FPGA
        clk_period = 1/frequency;
        
        % Number of bits for the input (x)
        sim_n_bit_in = circuit_parameters.inputResolution;
            
        % Number of bits for the output (y)
        sim_n_bit_out = circuit_parameters.outputResolution;
        
        if circuit_parameters.useADC == 1
            sim_bin_point_in = 0;
        else
            sim_bin_point_in = sim_n_bit_in - ceil(log2(max(abs([range.umax, range.umin, range.pmax, range.pmin, range.ymax, range.ymin]),[],'all'))) - 1;
        end

        if circuit_parameters.useDAC == 1
            sim_bin_point_out = 0;
        else
            sim_bin_point_out = sim_n_bit_out - ceil(log2(max(abs([range.xmax, range.xmin, range.dmax, range.dmin]),[],'all'))) - 1;
        end        
        
        % Compute arrays to transform the outputs from their circuit range to the
        % actual (model) range: u = sim_u_scale_gain.*u + sim_u_scale_bias    
        sim_u_scale_gain = (circuit_parameters.inputRange.max(1)-circuit_parameters.inputRange.min(1))./...
            (range.umax(:)-range.umin(:));
        sim_u_scale_bias = (range.umax(:)+range.umin(:))/2-...
            ((circuit_parameters.inputRange.max(1)+circuit_parameters.inputRange.min(1))/2./sim_u_scale_gain);
        
        if np == 0
            sim_p_scale_gain = [];
            sim_p_scale_bias = [];
        else
            sim_p_scale_gain = (circuit_parameters.inputRange.max(1)-circuit_parameters.inputRange.min(1))./...
                (range.pmax(:)-range.pmin(:));
            sim_p_scale_bias = (range.pmax+range.pmin)/2-...
                ((circuit_parameters.inputRange.max(1)+circuit_parameters.inputRange.min(1))/2./sim_p_scale_gain);
        end
        
        % Compute arrays to transform the inputs from their actual (model) range to
        % the circuit range: y_cir = sim_y_scale_gain.*y + sim_y_scale_bias
        sim_y_scale_gain = (circuit_parameters.inputRange.max(1)-circuit_parameters.inputRange.min(1))./...
            (range.ymax(:)-range.ymin(:));
        sim_y_scale_bias = (range.ymax(:)+range.ymin(:))/2-...
            ((circuit_parameters.inputRange.max(1)+circuit_parameters.inputRange.min(1))/2./sim_y_scale_gain);
        
        % Compute arrays to transform the inputs from their actual (model) range to
        % the circuit range: x_cir = sim_x_scale_gain.*x + sim_x_scale_bias
        sim_x_scale_gain = (range.xmax(:)-range.xmin(:))./(circuit_parameters.outputRange.max(1)-circuit_parameters.outputRange.min(1));
        sim_x_scale_bias = (range.xmax(:)+range.xmin(:))/2-...
            ((circuit_parameters.outputRange.max(1)+circuit_parameters.outputRange.min(1))/2.*sim_x_scale_gain);        
        
        if nd == 0
            sim_d_scale_gain = [];
            sim_d_scale_bias = [];
        else
            sim_d_scale_gain = (range.dmax(:)-range.dmin(:))./(circuit_parameters.outputRange.max(1)-circuit_parameters.outputRange.min(1));
            sim_d_scale_bias = (range.dmax(:)+range.dmin(:))/2-...
            ((circuit_parameters.outputRange.max(1)+circuit_parameters.outputRange.min(1))/2.*sim_d_scale_gain);    
        end
        
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
            'sim_nx','sim_ny','sim_np','sim_nd','sim_u_scale_gain','sim_u_scale_bias',...
            'sim_x_scale_gain','sim_x_scale_bias','sim_p_scale_gain','sim_p_scale_bias',...
            'sim_d_scale_gain','sim_d_scale_bias','sim_y_scale_gain','sim_y_scale_bias',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames','clk_period',...
            'sim_n_bit_in','sim_n_bit_out','sim_nu','sim_bin_point_in','sim_bin_point_out');
        
        preFolder = ['circuit_model/',circuit_parameters.inputRepresentation,'_in/'];
        
        % Add System Generator token
        add_block('sysgenController/ System Generator', [model, '/ System Generator']);

        % Add black box
        add_block('sysgenController/Controller', [model, '/Observer']);

        if nu ~= 0
            % Add monitor for output (u)
            add_block('simulink/Sinks/Scope', [model, '/Monitor u']);
            set_param([model, '/Monitor u'], 'BackgroundColor', 'Cyan');
            set_param([model, '/Monitor u'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_U', 'DataLoggingSaveFormat', 'Structure With Time');
            add_line(model, 'Constant input/1', 'Monitor u/1', 'autorouting', 'smart');
        end

        % Add reset input (step)
        add_block('sysgenController/bool_input', [model, '/reset']);
        add_block('simulink/Sources/Step', [model, '/rst']);
        set_param([model, '/rst'], 'Time', 'clk_period/2', 'Before', '0', 'After', '1');
        add_line(model, 'rst/1', 'reset/1', 'autorouting', 'smart');
        add_line(model, 'reset/1', 'Observer/1', 'autorouting', 'smart');

        if circuit_parameters.useADC == 1
            
            % Scale input (u)
            add_block('simulink/Ports & Subsystems/Subsystem', [model, '/controlScale']);
            set_param([model, '/controlScale'], 'ContentPreviewEnabled', 'off');
            delete_line([model, '/controlScale'], 'In1/1', 'Out1/1');
            set_param([model, '/controlScale/In1'], 'Name', 'u');
            set_param([model, '/controlScale/Out1'], 'Name', 'u_cir');
            add_block('simulink/Math Operations/Add', [model, '/controlScale/Add']);
            set_param([model, '/controlScale/Add'], 'Inputs', '+-');
            add_block('simulink/Math Operations/Product', [model, '/controlScale/Product']);
            add_block('simulink/Sources/Constant', [model, '/controlScale/bias']);
            set_param([model, '/controlScale/bias'], 'Value', 'sim_u_scale_bias');
            add_block('simulink/Sources/Constant', [model, '/controlScale/gain']);
            set_param([model, '/controlScale/gain'], 'Value', 'sim_u_scale_gain');
            add_line([model, '/controlScale'], 'u/1', 'Add/1', 'autorouting', 'smart');
            add_line([model, '/controlScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
            add_line([model, '/controlScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
            add_line([model, '/controlScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
            add_line([model, '/controlScale'], 'Product/1', 'u_cir/1', 'autorouting', 'smart');
            Simulink.BlockDiagram.arrangeSystem([model, '/controlScale']);
            add_line(model, 'Constant input/1', 'controlScale/1', 'autorouting', 'smart');
        
            % Scale parameters (p)
            if np ~= 0
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
        
            % Scale input (y)
            add_block('simulink/Ports & Subsystems/Subsystem', [model, '/outputScale']);
            set_param([model, '/outputScale'], 'ContentPreviewEnabled', 'off');
            delete_line([model, '/outputScale'], 'In1/1', 'Out1/1');
            set_param([model, '/outputScale/In1'], 'Name', 'y');
            set_param([model, '/outputScale/Out1'], 'Name', 'y_cir');
            add_block('simulink/Math Operations/Add', [model, '/outputScale/Add']);
            set_param([model, '/outputScale/Add'], 'Inputs', '+-');
            add_block('simulink/Math Operations/Product', [model, '/outputScale/Product']);
            add_block('simulink/Sources/Constant', [model, '/outputScale/bias']);
            set_param([model, '/outputScale/bias'], 'Value', 'sim_y_scale_bias');
            add_block('simulink/Sources/Constant', [model, '/outputScale/gain']);
            set_param([model, '/outputScale/gain'], 'Value', 'sim_y_scale_gain');
            add_line([model, '/outputScale'], 'y/1', 'Add/1', 'autorouting', 'smart');
            add_line([model, '/outputScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
            add_line([model, '/outputScale'], 'Add/1', 'Product/1', 'autorouting', 'smart');
            add_line([model, '/outputScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
            add_line([model, '/outputScale'], 'Product/1', 'y_cir/1', 'autorouting', 'smart');
            Simulink.BlockDiagram.arrangeSystem([model, '/outputScale']);
            add_line(model, 'LTI System/2', 'outputScale/1', 'autorouting', 'smart');

            % Scale output (x)
            add_block('simulink/Ports & Subsystems/Subsystem', [model, '/stateScale']);
            set_param([model, '/stateScale'], 'ContentPreviewEnabled', 'off');
            delete_line([model, '/stateScale'], 'In1/1', 'Out1/1');
            set_param([model, '/stateScale/In1'], 'Name', 'x_cir');
            set_param([model, '/stateScale/Out1'], 'Name', 'x');
            add_block('simulink/Math Operations/Product', [model, '/stateScale/Product']);
            add_block('simulink/Math Operations/Add', [model, '/stateScale/Add']);
            add_block('simulink/Sources/Constant', [model, '/stateScale/gain']);
            set_param([model, '/stateScale/gain'], 'Value', 'sim_x_scale_gain');
            add_block('simulink/Sources/Constant', [model, '/stateScale/bias']);
            set_param([model, '/stateScale/bias'], 'Value', 'sim_x_scale_bias');
            add_line([model, '/stateScale'], 'x_cir/1', 'Product/1', 'autorouting', 'smart');
            add_line([model, '/stateScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
            add_line([model, '/stateScale'], 'Product/1', 'Add/1', 'autorouting', 'smart');
            add_line([model, '/stateScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
            add_line([model, '/stateScale'], 'Add/1', 'x/1', 'autorouting', 'smart');
            Simulink.BlockDiagram.arrangeSystem([model, '/stateScale']);

            % Scale unmeasurable inputs or parameters (d)
            if nd ~= 0  
                add_block('simulink/Ports & Subsystems/Subsystem', [model, '/unmeasInputScale']);
                set_param([model, '/unmeasInputScale'], 'ContentPreviewEnabled', 'off');
                delete_line([model, '/unmeasInputScale'], 'In1/1', 'Out1/1');
                set_param([model, '/unmeasInputScale/In1'], 'Name', 'd_cir');
                set_param([model, '/unmeasInputScale/Out1'], 'Name', 'd');
                add_block('simulink/Math Operations/Product', [model, '/unmeasInputScale/Product']);
                add_block('simulink/Math Operations/Add', [model, '/unmeasInputScale/Add']);
                add_block('simulink/Sources/Constant', [model, '/unmeasInputScale/gain']);
                set_param([model, '/unmeasInputScale/gain'], 'Value', 'sim_d_scale_gain');
                add_block('simulink/Sources/Constant', [model, '/unmeasInputScale/bias']);
                set_param([model, '/unmeasInputScale/bias'], 'Value', 'sim_d_scale_bias');
                add_line([model, '/unmeasInputScale'], 'd_cir/1', 'Product/1', 'autorouting', 'smart');
                add_line([model, '/unmeasInputScale'], 'gain/1', 'Product/2', 'autorouting', 'smart');
                add_line([model, '/unmeasInputScale'], 'Product/1', 'Add/1', 'autorouting', 'smart');
                add_line([model, '/unmeasInputScale'], 'bias/1', 'Add/2', 'autorouting', 'smart');
                add_line([model, '/unmeasInputScale'], 'Add/1', 'd/1', 'autorouting', 'smart');
                Simulink.BlockDiagram.arrangeSystem([model, '/unmeasInputScale']);
            end
            
        end
        
        % Define inputs and outputs of the observer
        if nu == 1
            if strcmp(circuit_parameters.inputRepresentation,'signed')
                add_block('sysgenController/signed_fixed_input', [model, '/u1']);
            else
                add_block('sysgenController/unsigned_fixed_input', [model, '/u1']);
            end
            if circuit_parameters.useADC == 1
                add_line(model, 'controlScale/1', 'u1/1', 'autorouting', 'smart');
            else
                add_line(model, '/1', 'u1/1', 'autorouting', 'smart');
            end
            add_line(model, 'u1/1', 'Observer/2', 'autorouting', 'smart');
        else
            add_block('simulink/Signal Routing/Demux', [model, '/DemuxU']);
            set_param([model, '/DemuxU'], 'Outputs', 'sim_nu');
            if circuit_parameters.useADC == 1
                add_line(model, 'controlScale/1', 'DemuxU/1', 'autorouting', 'smart');
            else
                add_line(model, 'Constant input/1', 'DemuxU/1', 'autorouting', 'smart');
            end
            for i=1:nu
                if strcmp(circuit_parameters.inputRepresentation,'signed')
                    add_block('sysgenController/signed_fixed_input', [model, '/u', num2str(i)]);
                else
                    add_block('sysgenController/unsigned_fixed_input', [model, '/u', num2str(i)]);
                end
                add_line(model, ['DemuxU/', num2str(i)], ['u', num2str(i), '/1'], 'autorouting', 'smart');
                add_line(model, ['u', num2str(i), '/1'], ['Observer/', num2str(i+1)], 'autorouting', 'smart');
            end
        end
        
        if np ~= 0
            if np == 1
                if strcmp(circuit_parameters.inputRepresentation,'signed')
                    add_block('sysgenController/signed_fixed_input', [model, '/p1']);
                else
                    add_block('sysgenController/unsigned_fixed_input', [model, '/p1']);
                end
                if circuit_parameters.useADC == 1
                    add_line(model, 'parScale/1', 'p1/1', 'autorouting', 'smart');
                else
                    add_line(model, 'Parameter switch/1', 'p1/1', 'autorouting', 'smart');
                end
                add_line(model, 'p1/1', ['Observer/', num2str(nu+2)], 'autorouting', 'smart');
            else
                add_block('simulink/Signal Routing/Demux', [model, '/DemuxP']);
                set_param([model, '/DemuxP'], 'Outputs', 'sim_np');
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
                    add_line(model, ['p', num2str(i), '/1'], ['Observer/', num2str(nu+i+1)], 'autorouting', 'smart');
                end
            end
        end
        
        if ny == 1
            if strcmp(circuit_parameters.inputRepresentation,'signed')
                add_block('sysgenController/signed_fixed_input', [model, '/y1']);
            else
                add_block('sysgenController/unsigned_fixed_input', [model, '/y1']);
            end
            if circuit_parameters.useADC == 1
                add_line(model, 'outputScale/1', 'y1/1', 'autorouting', 'smart');
            else
                add_line(model, 'LTI System/2', 'y1/1', 'autorouting', 'smart');
            end
            add_line(model, 'y1/1', ['Observer/', num2str(nu+np+2)], 'autorouting', 'smart');
        else
            add_block('simulink/Signal Routing/Demux', [model, '/DemuxY']);
            set_param([model, '/DemuxY'], 'Outputs', 'sim_ny');
            if circuit_parameters.useADC == 1
                add_line(model, 'outputScale/1', 'DemuxY/1', 'autorouting', 'smart');
            else
                add_line(model, 'LTI System/2', 'DemuxY/1', 'autorouting', 'smart');
            end
            for i=1:ny
                if strcmp(circuit_parameters.inputRepresentation,'signed')
                    add_block('sysgenController/signed_fixed_input', [model, '/y', num2str(i)]);
                else
                    add_block('sysgenController/unsigned_fixed_input', [model, '/y', num2str(i)]);
                end
                add_line(model, ['DemuxY/', num2str(i)], ['y', num2str(i), '/1'], 'autorouting', 'smart');
                add_line(model, ['y', num2str(i), '/1'], ['Observer/', num2str(nu+np+i+1)], 'autorouting', 'smart');
            end
        end
        
        if nx == 1
            add_block('xbsBasic_r4/Gateway Out', [model, '/x1_stim']);
            add_line(model, 'Observer/1', 'x1_stim/1', 'autorouting', 'smart');
            if circuit_parameters.useDAC == 1
                add_line(model, 'x1_stim/1', 'stateScale/1', 'autorouting', 'smart');
                add_line(model, 'stateScale/1', 'Monitor xest/1', 'autorouting', 'smart');
            else
                add_line(model, 'x1_stim/1', 'Monitor xest/1', 'autorouting', 'smart');
            end
        else
            add_block('simulink/Signal Routing/Mux', [model, '/MuxX']);
            set_param([model, '/MuxX'], 'Inputs', num2str(nx));
            if circuit_parameters.useDAC == 1
                add_line(model, 'MuxX/1', 'stateScale/1', 'autorouting', 'smart');
                add_line(model, 'stateScale/1', 'Monitor xest/1', 'autorouting', 'smart');
            else
                add_line(model, 'MuxX/1', 'Monitor xest/1', 'autorouting', 'smart');
            end
            for i=1:nx
                add_block('xbsBasic_r4/Gateway Out', [model, '/x', num2str(i), '_stim']);
                add_line(model, ['Observer/', num2str(i)], ['x', num2str(i), '_stim/1'], 'autorouting', 'smart');
                add_line(model, ['x', num2str(i), '_stim/1'], ['MuxX/', num2str(i)], 'autorouting', 'smart');
            end
        end

        if nd ~= 0
            add_block('simulink/Sinks/Scope', [model, '/Monitor dest']);
            set_param([model, '/Monitor dest'], 'BackgroundColor', 'Cyan');
            set_param([model, '/Monitor dest'], 'DataLoggingDecimateData', 'On','DataLoggingDecimation', '1', 'DataLogging', 'On', 'DataLoggingVariableName', 'res_Dest', 'DataLoggingSaveFormat', 'Structure With Time');
            if nd == 1
                add_block('xbsBasic_r4/Gateway Out', [model, '/d1_stim']);
                add_line(model, ['Observer/', num2str(nx+1)], 'd1_stim/1', 'autorouting', 'smart');
                if circuit_parameters.useDAC == 1
                    add_line(model, 'd1_stim/1', 'unmeasInputScale/1', 'autorouting', 'smart');
                    add_line(model, 'unmeasInputScale/1', 'Monitor dest/1', 'autorouting', 'smart');
                else
                    add_line(model, 'd1_stim/1', 'Monitor dest/1', 'autorouting', 'smart');
                end
            else
                add_block('simulink/Signal Routing/MuxD', [model, '/MuxD']);
                set_param([model, '/MuxD'], 'Inputs', num2str(nd));
                if circuit_parameters.useDAC == 1
                    add_line(model, 'MuxD/1', 'unmeasInputScale/1', 'autorouting', 'smart');
                    add_line(model, 'unmeasInputScale/1', 'Monitor dest/1', 'autorouting', 'smart');
                else
                    add_line(model, 'MuxD/1', 'Monitor dest/1', 'autorouting', 'smart');
                end
                for i=1:nd
                    add_block('xbsBasic_r4/Gateway Out', [model, '/d', num2str(i), '_stim']);
                    add_line(model, ['Observer/', num2str(i + nx)], ['d', num2str(i), '_stim/1'], 'autorouting', 'smart');
                    add_line(model, ['x', num2str(i), '_stim/1'], ['MuxD/', num2str(i)], 'autorouting', 'smart');
                end
            end
        end
        
        % Add outputs 'sample'
        add_block('xbsBasic_r4/Gateway Out', [model, '/sample']);
        add_block('simulink/Sinks/Scope', [model, '/Monitor sample']);
        add_line(model, ['Observer/', num2str(nx+nd+1)], 'sample/1', 'autorouting', 'smart');
        add_line(model, 'sample/1', 'Monitor sample/1', 'autorouting', 'smart');
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
    
    if ~simulateVHDL
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
    else
        fprintf(fp,'%% Observer initial condition for the unmeasurable inputs and states\n %% are set in VHDL\n');
    end
    
    fprintf(fp,'%% Simulation time\n');
    fprintf(fp,'sim_T = 1;\n');
    fprintf(fp,'\n');
    fprintf(fp,'%% If you do not use the signal builder for the inputs, parameters and \n');
    fprintf(fp,'%% unmeasurable inputs, you can set their constant values below\n');
    fprintf(fp,'\n');
    
    if nu > 0
        fprintf(fp,'%% Constant input\n');
        fprintf(fp,'sim_u = [');
        for i = 1:object.nu-1
            fprintf(fp,'0;');
        end
        fprintf(fp,'0];\n');
        fprintf(fp,'\n');
    end
    
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

