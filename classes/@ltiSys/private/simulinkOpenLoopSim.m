function simulinkOpenLoopSim(object,options)
% simulinkOpenLoopSim          Generates the Simulink model of the 
%                              open loop dynamical system

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

if strcmpi(choice,'Overwrite')    
    
    % Create folder
    mkdir(options.folder)
    
    cd(options.folder)
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    sim_np = np;
    sim_nu = nu;
   
    % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D',...
        'sim_nx','sim_ny','sim_np','sim_nd','sim_xnames','sim_unames','sim_pnames',...
        'sim_dnames','sim_ynames');
    
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
            add_line([model, '/LTI System'], 'Const/1', 'Mux/2', 'autorouting', 'smart');
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

