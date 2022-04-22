%% PLOT SCOPE
% Plots the results of Simulink simulations

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------
   
% Number of states
nx = size(res_X.signals.values,2);

if exist('sim_ctrl','var')
    tracking = sim_ctrl.isTracking();
    trackvar = sim_ctrl.getTrackingVariable();
else
    tracking = 0;
end

% Plot states
figure('Name','System states')
for i = 1:nx
    subplot(nx,1,i)
    hold on
    plot(res_X.time,res_X.signals.values(:,i),'b')
    if exist('res_Xest','var')
        stairs(res_Xest.time,res_Xest.signals.values(:,i),'r')
        legend('real','estimated')
    end
    if tracking && any(trackvar == i)
        plot(res_Xref.time,res_Xref.signals.values,'k')
    end
    xlabel('t')
    ylabel(sim_xnames{i})
    grid on
end
    
% Plot inputs
if exist('res_U','var')
    nu = size(res_U.signals.values,2);
    figure('Name','System inputs')
    for i = 1:nu
        subplot(nu,1,i)
        hold on
        stairs(res_U.time,res_U.signals.values(:,i),'b')
        xlabel('t')
        ylabel(sim_unames{i})
        grid on
    end
end

% Plot outputs
if exist('res_Y','var')
    ny = size(res_Y.signals.values,2);
    figure('Name','System outputs')
    for i = 1:ny
        subplot(ny,1,i)
        hold on
        plot(res_Y.time,res_Y.signals.values(:,i),'b')
        if exist('res_Yest','var')
            stairs(res_Yest.time,res_Yest.signals.values(:,i),'r')
            legend('real','estimated')
        end
        xlabel('t')
        ylabel(sim_ynames{i})
        grid on
    end
end

% Plot parameters
if exist('res_P','var')
    np = size(res_P.signals.values,2);
    figure('Name','System parameters')
    for i = 1:np
        subplot(np,1,i)
        hold on
        plot(res_P.time,res_P.signals.values(:,i),'b')
        xlabel('t')
        ylabel(sim_pnames{i})
        grid on
    end
end

% Plot unmeasurable inputs
if exist('res_D','var')
    nd = size(res_D.signals.values,2);
    figure('Name','System unmeasurable inputs')
    for i = 1:nd
        subplot(nd,1,i)
        hold on
        plot(res_D.time,res_D.signals.values(:,i),'b')
        if exist('res_Dest','var')
            stairs(res_Dest.time,res_Dest.signals.values(:,i),'r')
            legend('real','estimated')
        end
        xlabel('t')
        ylabel(sim_dnames{i})
        grid on
    end
end