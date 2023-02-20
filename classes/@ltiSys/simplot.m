function signals = simplot(object,T,x0,varargin)
% simplot    Simulates the dynamical system and plots results
%
% Performs a simulation by means of method ltiSys/sim and plots the
% results. The syntax is the same of method ltiSys/sim.
% SIGNALS = simplot(OBJ,T,X0)
% SIGNALS = simplot(OBJ,T,X0,OPTS)
% SIGNALS = simplot(OBJ,T,X0,REF)
% SIGNALS = simplot(OBJ,T,X0,REF,OPTS)
% SIGNALS = simplot(OBJ,T,X0,P,D)
% SIGNALS = simplot(OBJ,T,X0,P,D,OPTS)
% SIGNALS = simplot(OBJ,T,X0,P,D,REF)
% SIGNALS = simplot(OBJ,T,X0,P,D,REF,OPTS)
%
% If a controller is associated to the system, and the controller contains
% the information about the constraintsis provided, the constraints are plotted in
% order to visually check if they are fulfilled. Hard constraints are
% plotted with solid black lines, soft constraints with dashed black lines.
%
% See also: ltiSys/sim.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

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

if isempty(varargin)
    signals = object.sim(T,x0);
else
    signals = object.sim(T,x0,varargin{:});
end

ctrl = object.getController();
if ~isempty(ctrl)
    info = ctrl.getInformation();
    constr = info.constr;
else
    constr = [];
end

% Retrieve names
xnames = object.getStateNames();
unames = object.getInputNames();
ynames = object.getOutputNames();
dnames = object.getUnmeasurableInputNames();


% Number of states, inputs and outputs
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;

nxref = numel(object.getController.getTrackingVariable);

npts = numel(signals.time);

if ~isempty(signals.time)
    
    if object.isContinuousTime
        
        % Extract constraints matrices at following following time instant
        if ~isempty(constr)
            N = constr.getTimeHorizon();
            if N >= 1
                [H, K] = constr.getAllConstraints('hard',1);
                [Hs, Ks] = constr.getAllConstraints('soft',1);
%                 if nxref ~= 0
%                     H = [H zeros(size(H,1),length(signals.ref))];
%                     K = [K zeros(size(K,1),length(signals.ref))];
%                 end
                
            else
                H = [];
                K = [];
                Hs = [];
                Ks = [];
            end
        else
            H = [];
            K = [];
            Hs = [];
            Ks = [];
        end
        
        
        % Plot states
        figure('Name','System states')
        for i = 1:nx
            subplot(nx,1,i)
            hold on
            
            if ~isempty(H)
                Hleft = H(:,i);
                Hright = H(:,setdiff(1:size(H,2),i));
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                if ~isempty(Hleft)
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    
                    if ~isempty(Krightp)
                        if nxref ~= 0
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        else
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        end
                        
                    end
                    if ~isempty(Krightn)
                        if nxref ~= 0
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')
                        else
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')
                        end
                        
                    end
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    
                    if ~isempty(Ksrightp)
                        if nxref ~= 0
                            up = (repmat(Ksrightp',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                                repmat(Hsleftp',npts,1);
                            plot(signals.time,min(up,[],2),'--k')
                        else
                            up = (repmat(Ksrightp',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hsrightp')./...
                                repmat(Hsleftp',npts,1);
                            plot(signals.time,min(up,[],2),'--k')
                        end
                    end
                    if ~isempty(Ksrightn)
                        if nxref ~= 0
                            down = (repmat(Ksrightn',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                                repmat(Hsleftn',npts,1);
                            plot(signals.time,max(down,[],2),'--k')
                        else
                            down = (repmat(Ksrightn',npts,1)-...
                                [signals.state(:,setdiff(1:nx,i))...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hsrightn')./...
                                repmat(Hsleftn',npts,1);
                            plot(signals.time,max(down,[],2),'--k')
                        end
                    end
                end
            end
            
            plot(signals.time,signals.state(:,i),'b')
            if object.hasObserver()
                plot(signals.time,signals.est_state(:,i),'r')
                legend('real','estimated')
            end
            if object.hasController
                ctrl = object.getController();
                if ctrl.isTracking
                    if strcmpi(ctrl.getControlVariable,'state')
                        trackvar = ctrl.getTrackingVariable();
                        idx = find(trackvar == i);
                        if ~isempty(idx) && signals.ref(idx,end) ~= signals.ref(idx,end)
                            plot([signals.time(1) signals.time(end)],[signals.ref(idx) signals.ref(idx)],'c')
                        end
                    end
                end
            end
            xlabel('t')
            ylabel(xnames{i})
            grid on
            
            % Slightly enlarge axes
            ax = axis;
            delta = 0.05*(ax(4)-ax(3));
            ax(3) = ax(3)-delta;
            ax(4) = ax(4)+delta;
            axis(ax);
            
        end
        
        % Plot unmeasurable inputs
        if object.hasObserver()
            if object.nd > 0
                figure('Name','Unmeasurable inputs')
                for i = 1:nd
                    subplot(nd,1,i)
                    hold on
                    plot([signals.time(1) signals.time(end)],...
                        [signals.unmeasurable_input(i) signals.unmeasurable_input(i)],'b')
                    plot(signals.time,signals.est_state(:,nx+i),'r')
                    xlabel('t')
                    ylabel(dnames{i})
                    grid on
                end
            end
        end
        
        % Extract constraints matrices at current time instant
        if ~isempty(constr)
            [H, K] = constr.getAllConstraints('hard',0);
            [Hs, Ks] = constr.getAllConstraints('soft',0);
%             if nxref ~= 0
%                 H = [H zeros(size(H,1),length(signals.ref))];
%                 K = [K zeros(size(K,1),length(signals.ref))];
%             end
        else
            H = [];
            Hs = [];
            K = [];
            Ks = [];
        end
        
        % Plot inputs
        figure('Name','System inputs')
        for i = 1:nu
            subplot(nu,1,i)
            hold on
            
            if ~isempty(H)
                
                Hleft = H(:,nx+i);
                Hright = H(:,setdiff(1:size(H,2),nx+i));
                
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                
                if ~isempty(Hleft)
                    
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    
                    if ~isempty(Krightp)
                        if nxref ~= 0
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state...
                                signals.input(:,setdiff(1:nu,i))...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        else
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state...
                                signals.input(:,setdiff(1:nu,i))...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        end
                    end
                    if ~isempty(Krightn)
                        if nxref ~= 0
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state...
                                signals.input(:,setdiff(1:nu,i))...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')
                        else
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state...
                                signals.input(:,setdiff(1:nu,i))...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')
                        end
                    end                    
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,nx+i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),nx+i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    
                    if ~isempty(Ksrightp)
                        up = (repmat(Ksrightp',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                            repmat(Hsleftp',npts,1);
                        plot(signals.time,min(up,[],2),'--k')
                    end
                    if ~isempty(Ksrightn)
                        down = (repmat(Ksrightn',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                            repmat(Hsleftn',npts,1);
                        plot(signals.time,max(down,[],2),'--k')
                    end
                    
                    
                end
            end
            
            plot(signals.time,signals.input(:,i),'b')
            xlabel('t')
            ylabel(unames{i})
            grid on
            
            % Slightly enlarge axes
            ax = axis;
            delta = 0.05*(ax(4)-ax(3));
            ax(3) = ax(3)-delta;
            ax(4) = ax(4)+delta;
            axis(ax);
        end
        
        % Plot outputs
        figure('Name','System outputs')
        for i = 1:ny
            subplot(ny,1,i)
            hold on

            if ~isempty(H)
                Hleft = H(:,nx+nu+np+nd+i);
                Hright = H(:,setdiff(1:size(H,2),nx+nu+np+nd+i));
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                if ~isempty(Hleft)
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    
                    if ~isempty(Krightp)
                        if nxref ~= 0
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output(:,setdiff(1:ny,i))...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        else
                            up = (repmat(Krightp',npts,1)-...
                                [signals.state...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output(:,setdiff(1:ny,i))]*Hrightp')./...
                                repmat(Hleftp',npts,1);
                            plot(signals.time,min(up,[],2),'k')
                        end
                    end
                    if ~isempty(Krightn)
                        if nxref ~= 0
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output(:,setdiff(1:ny,i))...
                                repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')
                        else
                            down = (repmat(Krightn',npts,1)-...
                                [signals.state...
                                signals.input...
                                repmat(signals.parameter',npts,1)...
                                repmat(signals.unmeasurable_input',npts,1)...
                                signals.output(:,setdiff(1:ny,i))]*Hrightn')./...
                                repmat(Hleftn',npts,1);
                            plot(signals.time,max(down,[],2),'k')                            
                        end
                    end
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,nx+nu+np+nd+i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),nx+nu+np+nd+i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    
                    if ~isempty(Ksrightp)
                        up = (repmat(Ksrightp',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                            repmat(Hsleftp',npts,1);
                        plot(signals.time,min(up,[],2),'--k')
                    end
                    if ~isempty(Ksrightn)
                        down = (repmat(Ksrightn',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                            repmat(Hsleftn',npts,1);
                        plot(signals.time,max(down,[],2),'--k')
                    end
                end
            end

            plot(signals.time,signals.output(:,i),'b')
            if object.hasObserver()
                plot(signals.time,signals.est_output(:,i),'r')
                legend('real','estimated')
            end
            if object.hasController
                ctrl = object.getController();
                if ctrl.isTracking
                    if strcmpi(ctrl.getControlVariable,'output')
                        trackvar = ctrl.getTrackingVariable();
                        idx = find(trackvar == i);
                        if ~isempty(idx)
                            plot([signals.time(1) signals.time(end)],[signals.ref(idx) signals.ref(idx)],'c')
                        end
                    end
                end
            end
            xlabel('t')
            ylabel(ynames{i})
            grid on
        end
        
    else
        
        % Extract constraints matrices at following time instant
        if ~isempty(constr)
            [H, K] = constr.getAllConstraints('hard',1);
            [Hs, Ks] = constr.getAllConstraints('soft',1);
            if nxref ~= 0
                H = [H zeros(size(H,1),length(signals.ref))];
                K = [K zeros(size(K,1),length(signals.ref))];
            end
        else
            H = [];
            Hs = [];
            K = [];
            Ks = [];
        end
        
        % Plot states
        figure('Name','System states')
        for i = 1:nx
            subplot(nx,1,i)
            hold on
            
            if ~isempty(H)
                Hleft = H(:,i);
                Hright = H(:,setdiff(1:size(H,2),i));
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                if ~isempty(Hleft)
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    if ~isempty(Krightp)
                        up = (repmat(Krightp',npts,1)-...
                            [signals.state(:,setdiff(1:nx,i))...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                            repmat(Hleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'k')
                    end
                    if ~isempty(Krightn)
                        down = (repmat(Krightn',npts,1)-...
                            [signals.state(:,setdiff(1:nx,i))...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                            repmat(Hleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'k')
                    end
                    
                    
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    
                    if ~isempty(Ksrightp)
                        up = (repmat(Ksrightp',npts,1)-...
                            [signals.state(:,setdiff(1:nx,i))...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                            repmat(Hsleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'--k')
                    end
                    if ~isempty(Ksrightn)
                        down = (repmat(Ksrightn',npts,1)-...
                            [signals.state(:,setdiff(1:nx,i))...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                            repmat(Hsleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'--k')
                    end
                    
                    
                end
            end
            
            stairs(signals.time,signals.state(:,i),'b')
            if object.hasObserver()
                stairs(signals.time,signals.est_state(:,i),'r')
                legend('real','estimated')
            end
            xlabel('t')
            ylabel(xnames{i})
            grid on
            
            % Slightly enlarge axes
            ax = axis;
            delta = 0.05*(ax(4)-ax(3));
            ax(3) = ax(3)-delta;
            ax(4) = ax(4)+delta;
            axis(ax);
        end
        
        % Plot unmeasurable inputs
        if object.hasObserver()
            figure('Name','Unmeasurable inputs')
            for i = 1:nd
                subplot(nd,1,i)
                hold on
                plot([signals.time(1) signals.time(end)],...
                    [signals.unmeasurable_input(i) signals.unmeasurable_input(i)],'b')
                plot(signals.time,signals.est_state(:,nx+i),'r')
                xlabel('t')
                ylabel(dnames{i})
                grid on
            end
        end
        
        % Extract constraints matrices at current time instant
        if ~isempty(constr)
            [H, K] = constr.getAllConstraints('hard',0);
            [Hs, Ks] = constr.getAllConstraints('soft',0);
            if nxref ~= 0
                H = [H zeros(size(H,1),length(signals.ref))];
                K = [K zeros(size(K,1),length(signals.ref))];
            end
        else
            H = [];
            Hs = [];
            K = [];
            Ks = [];
        end
        
        % Plot inputs
        figure('Name','System inputs')
        for i = 1:nu
            subplot(nu,1,i)
            hold on
            
            if ~isempty(H)
                Hleft = H(:,nx+i);
                Hright = H(:,setdiff(1:size(H,2),nx+i));
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                if ~isempty(Hleft)
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    if ~isempty(Krightp)
                        up = (repmat(Krightp',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                            repmat(Hleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'k')
                    end
                    if ~isempty(Krightn)
                        down = (repmat(Krightn',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                            repmat(Hleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'k')
                    end
                    
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,nx+i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),nx+i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    if ~isempty(Ksrightp)
                        up = (repmat(Ksrightp',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                            repmat(Hsleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'--k')
                    end
                    if ~isempty(Ksrightn)
                        down = (repmat(Ksrightn',npts,1)-...
                            [signals.state...
                            signals.input(:,setdiff(1:nu,i))...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                            repmat(Hsleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'--k')
                    end
                    
                    
                end
            end
            
            stairs(signals.time,signals.input(:,i),'b')
            xlabel('t')
            ylabel(unames{i})
            grid on
            
            % Slightly enlarge axes
            ax = axis;
            delta = 0.05*(ax(4)-ax(3));
            ax(3) = ax(3)-delta;
            ax(4) = ax(4)+delta;
            axis(ax);
            
        end
        
        % Plot outputs
        figure('Name','System outputs')
        for i = 1:ny
            subplot(ny,1,i)
            hold on

            if ~isempty(H)
                Hleft = H(:,nx+nu+np+nd+i);
                Hright = H(:,setdiff(1:size(H,2),i));
                idx = Hleft~=0;
                Hleft = Hleft(idx);
                if ~isempty(Hleft)
                    Hright = Hright(idx,:);
                    Kright = K(idx);
                    
                    % Separate positive and negative constraints
                    idxp = Hleft > 0;
                    idxn = Hleft < 0;
                    
                    Krightp = Kright(idxp);
                    Krightn = Kright(idxn);
                    Hrightp = Hright(idxp,:);
                    Hrightn = Hright(idxn,:);
                    Hleftp = Hleft(idxp);
                    Hleftn = Hleft(idxn);
                    if ~isempty(Krightp)
                        up = (repmat(Krightp',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightp')./...
                            repmat(Hleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'k')
                    end
                    if ~isempty(Krightn)
                        down = (repmat(Krightn',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hrightn')./...
                            repmat(Hleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'k')
                    end
                    
                    
                end
            end
            if ~isempty(Hs)
                Hsleft = Hs(:,nx+nu+np+nd+i);
                Hsright = Hs(:,setdiff(1:size(Hs,2),i));
                idx = Hsleft~=0;
                Hsleft = Hsleft(idx);
                if ~isempty(Hsleft)
                    Hsright = Hsright(idx,:);
                    Ksright = Ks(idx);
                    % Separate positive and negative constraints
                    idxp = Hsleft > 0;
                    idxn = Hsleft < 0;
                    
                    Ksrightp = Ksright(idxp);
                    Ksrightn = Ksright(idxn);
                    Hsrightp = Hsright(idxp,:);
                    Hsrightn = Hsright(idxn,:);
                    Hsleftp = Hsleft(idxp);
                    Hsleftn = Hsleft(idxn);
                    
                    if ~isempty(Ksrightp)
                        up = (repmat(Ksrightp',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightp')./...
                            repmat(Hsleftp',npts,1);
                        stairs(signals.time,min(up,[],2),'--k')
                    end
                    if ~isempty(Ksrightn)
                        down = (repmat(Ksrightn',npts,1)-...
                            [signals.state...
                            signals.input...
                            repmat(signals.parameter',npts,1)...
                            repmat(signals.unmeasurable_input',npts,1)...
                            signals.output(:,setdiff(1:ny,i))...
                            repmat(signals.ref(:,end)',npts,1)]*Hsrightn')./...
                            repmat(Hsleftn',npts,1);
                        stairs(signals.time,max(down,[],2),'--k')
                    end
                    
                    
                end
            end

            stairs(signals.time,signals.output(:,i),'b')
            if object.hasObserver()
                stairs(signals.time,signals.est_output(:,i),'r')
                legend('real','estimated')
            end
            xlabel('t')
            ylabel(ynames{i})
            grid on
        end
        
    end
    
%     if object.hasController
    if ~isa(object.getController,'implicitMPCctrl')
        
        ctrl = object.getController();
        ndim = ctrl.getNumberOfDimensions;
        
        % Domain dimensions of the controller
        
        if ndim == 2
            object.controller.plotPartition();
            hold on
            plot(signals.state(:,1),signals.state(:,2),'k','linewidth',2)
            if object.hasObserver()
                plot(signals.est_state(:,1),signals.est_state(:,2),'--k','linewidth',2)
            end
        end
    end
    
end

