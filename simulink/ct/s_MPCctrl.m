function [sys,x0,str,ts] = s_MPCctrl(t,x,u,flag,ctrl)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts] = mdlInitializeSizes(ctrl);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
%     sys=mdlDerivatives(t,x,u,A,B,C,D);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys = mdlOutputs(t,u,ctrl);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(ctrl)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = ctrl.getNumberOfInputs();
sizes.NumInputs      = ctrl.getNumberOfDimensions();
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [ctrl.getSamplingTime() 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
% function sys=mdlDerivatives(t,x,u,A,B,C,D)
% 
% sys = A*x + B*u;

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys = mdlOutputs(t,u,ctrl)

nx = ctrl.getNumberOfStates();
np = ctrl.getNumberOfParameters();
nd = ctrl.getNumberOfUnmeasurableInputs();

% x = u(1:nx)';
% p = u(nx+1:nx+np)';
% d = u(nx+np+1:nx+np+nd)';
x = u(1:nx);
p = u(nx+1:nx+np);
d = u(nx+np+1:nx+np+nd);

if ctrl.isTracking
    xref = u(nx+np+nd+1:end);
    sys = ctrl.eval(x,p,d,xref');
else
    sys = ctrl.eval(x,p,d);
end

% Backup controller
% If 0, when MPC controller returns NaN, an error is thrown
% If 1, when MPC controller returns NaN, the control function is set to 0
backup = 0;

if isnan(sys)
    if ~backup
        if t > 0
            error(['MPC controller returned NaN. System state may be outside ',...
                'the feasibility domain. Open file s_MPCctrl.m at line 94',...
                ' to enable a backup controller'])
        else
            % For t = 0, sometimes it happens that the input u is wrong
            % (maybe for some transient of the integration algorithms) and
            % the controller can be therefore infeasible. In this case it
            % is set to 0 without throwing errors.
            sys = zeros(size(sys));
        end
    else
        warning(['t = ',num2str(t),': backup controller (u = 0)']);
        sys = zeros(size(sys));
    end
end
   
% end mdlOutputs
