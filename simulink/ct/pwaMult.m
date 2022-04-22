function [sys,x0,str,ts]=pwaMult(t,x,u,flag,nx,np,nd,nu,H,K,C,D)

switch flag,
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes(np,nu,C,D);

        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3
        sys = mdlOutputs(u,nx,nu,np,nd,H,K,C,D);
        
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
function [sys,x0,str,ts] = mdlInitializeSizes(np,nu,C,D)
sizes = simsizes;
sizes.NumContStates = 0;
sizes.NumDiscStates = 0;
sizes.NumOutputs = size(C{1},1);
sizes.NumInputs = size([C{1} D{1}(:,nu+1:end)],2);

DirFeedthrough = 1;

sizes.DirFeedthrough = DirFeedthrough;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0 = [];
str = [];
ts = [];
% end mdlInitializeSizes
%



%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(u,nx,nu,np,nd,H,K,C,D)
for i=1:numel(C)    
    if all(H{i}*[u(1:nx);u(nx+nd+1:nx+nd+np);u(nx+1:nx+nd)] <= K{i})
        sys = [C{i} D{i}]*[u(1:nx+nd) ; zeros(nu,1); u(nx+nd+1:end)];
        break;
    end;
end;
% end mdlOutputs