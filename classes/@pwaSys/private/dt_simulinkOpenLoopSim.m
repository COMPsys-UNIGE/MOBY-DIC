function dt_simulinkOpenLoopSim(object,options)
% dt_simulinkOpenLoopSim    Generates a Simulink model for the 
%                           simulation of a discrete-time open-loop 
%                           dynamical system
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


% Default value
choice = 'Overwrite';

% If the model already existe, ask what to do...
if exist([options.folder,'pwaSys_model.mdl'],'file')
    
    % Construct a questdlg with three options
    choice = questdlg(['A model already exists in ',options.folder,'. Do you want to overwrite it completely or just update the system matrices?'], ...
        ' ', ...
        'Overwrite','Update','Cancel','Cancel');
    
end

% Extract matrices
[A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy,sim_H,sim_K] = object.getMatrices();

% Sample time
sim_Ts = object.Ts;

% Number of variables
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;

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

if strcmpi(choice,'Overwrite')
    
    % Create folder
    mkdir(options.folder)
    
    % Choose model to generate
    if nu > 0
        if np > 0
            if nd > 0
                model = 'pwaSys_model_upd.mdl';
            else
                model = 'pwaSys_model_up.mdl';
            end
        else
            if nd > 0
                model = 'pwaSys_model_ud.mdl';
            else
                model = 'pwaSys_model_u.mdl';
            end
        end
    else
        if np > 0
            if nd > 0
                model = 'pwaSys_model_pd.mdl';
            else
                model = 'pwaSys_model_p.mdl';
            end
        else
            if nd > 0
                model = 'pwaSys_model_d.mdl';
            else
                model = 'pwaSys_model.mdl';
            end
        end
    end
    
    % Copy model templates
    copyfile([options.simulinkFolder,'dt/',model],[options.folder,'pwaSys_model.mdl']);
    copyfile([options.simulinkFolder,'dt/','plotScope.m'],[options.folder,'plotScope.m']);
      copyfile([options.simulinkFolder,'dt/','s_pwaSys_d.m'],[options.folder,'s_pwaSys_d.m']);
   
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    
    % Save .mat files with the matrices
    save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D','sim_Ts',...
        'sim_H','sim_K','sim_nx','sim_ny','sim_nd',...
        'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames');
    
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
        fprintf(fp,'%% Constant inputs\n');
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
    
    sim_nx = nx;
    sim_ny = ny;
    sim_nd = nd;
    
    try
        % Save .mat files with the matrices
        save([options.folder,'simulink.mat'],'sim_A','sim_B','sim_C','sim_D','sim_Ts',...
            'sim_H','sim_K','sim_nx','sim_ny','sim_nd',...
            'sim_xnames','sim_unames','sim_pnames','sim_dnames','sim_ynames')
        
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






