function [sys,x0,str,ts]=s_pwaSys_d(t,x,u,flag,A,B,C,D,H,K,sim_x0,Ts)

switch flag,
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes(A,B,C,D,H,K,sim_x0,Ts);
        
        %%%%%%%%%%%%%%%
        % Update      %
        %%%%%%%%%%%%%%%
    case 2
        sys = mdlUpdate(t,x,u,A,B,C,D,H,K);
        
        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3
        sys = mdlOutputs(t,x,u,A,B,C,D,H,K);
        
        %%%%%%%%%%%%%
        % Terminate %
        %%%%%%%%%%%%%
    case 9
        sys = []; % do nothing
end
%end simom2
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes(A,B,C,D,H,K,sim_x0,Ts)
sizes = simsizes;
sizes.NumContStates = 0;
sizes.NumDiscStates = size(A{1},1);
sizes.NumOutputs = size(C{1},1);
sizes.NumInputs = size(B{1},2);

DirFeedthrough = 1;

nDyn = numel(A);

% for i=1:nDyn
%     if any(any(D{i}~=0))
%         DirFeedthrough = 1;
%     end
% end

sizes.DirFeedthrough = DirFeedthrough;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0 = sim_x0;
str = [];
ts = Ts;
% end mdlInitializeSizes
%

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys = mdlUpdate(t,x,u,A,B,C,D,H,K)
nx = size(A{1},1);
np_plus_nd = size(H{1},2)-nx;
uu = u(end-np_plus_nd:end-1);
for i=1:numel(A)
    if all(H{i}*[x;uu] <= K{i})
        sys = A{i}*x + B{i}*u;
        break;
    end;
end;
% end mdlDerivatives
%

%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u,A,B,C,D,H,K)
nx = size(A{1},1);
np_plus_nd = size(H{1},2)-nx;
uu = u(end-np_plus_nd:end-1);
for i=1:numel(A)
    if all(H{i}*[x;uu] <= K{i})
            sys = C{i}*x + D{i}*u;
        break;
    end;
end;
% end mdlOutputs