function object = setController(object,controller)
% setController   Associates a controller to the system
%
% OBJ = setController(OBJ,CTRL)
% CTRL must be a controller object. OBJ is the dynSys object.
%
% The observer must have the same number of inputs, states, parameters,
% and unmeasurable inputs. If the system is discrete-time, the
% sampling time of the controller must be equal or multiple of the system sampling
% time. If an observer is already associated to the system, the sampling
% time of the controller must be equal or multiple of the observer sampling time.
%
% See also dynSys/setController.
 
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

if object.nu == 0
    error('System has no inputs! Impossible to associate a controller.');
end

if ~isa(controller,'controller')
    error('CTRL must be a controller object')
end

% % TO DO: add automatic model augmenting when the controller is made for tracking 
% if isa(controller, 'explicitMPCctrl') && controller.getInformation.options.tracking
%     object = object.deltau();
% end

% Check compatibility
if controller.getNumberOfStates() ~= object.nx
    error('Controller and system have different number of states');
end
if controller.getNumberOfInputs() ~= object.nu
    error('Controller and system have different number of inputs');
end
if controller.getNumberOfParameters() ~= object.np
    error('Controller and system have different number of parameters');
end
if controller.getNumberOfUnmeasurableInputs() ~= object.nd
    error('Controller and system have different number of unmeasurable inputs');
end

% Retrieve rounding tolerance
pars = getMOBYDICpars();
tol = pars.roundtol;

% Check sampling time
if object.isDiscreteTime()
    if controller.getSamplingTime() < object.Ts
        error('Sampling time of the controller must not be lower than system sampling time');
    end
    
    % Ratio between sampling times
    Tsratio = controller.getSamplingTime()/object.Ts;
    
    % Round sampling time ratio, to prevent numerical issues
    Tsratio = round(Tsratio/tol)*tol;
    
    if floor(Tsratio) ~= Tsratio
        error('The sampling time of the controller must be a multiple of the system sampling time');
    end
end

% Check sampling time
if object.hasObserver()
    obs = object.observer;
    if controller.getSamplingTime() < obs.getSamplingTime()
        error('Sampling time of the controller must not be lower than observer sampling time');
    end
    
    % Ratio between sampling times
    Tsratio = controller.getSamplingTime()/obs.getSamplingTime();
    
    % Round sampling time ratio, to prevent numerical issues
    Tsratio = round(Tsratio/tol)*tol;
    
    if floor(Tsratio) ~= Tsratio
        error('The sampling time of the controller must be a multiple of the observer sampling time');
    end
    
    if isa(obs,'kalmanFilter') && isa(object,'ltiSys')
        D = object.getMatrices('D');
        if any(D ~= 0)
            disp('WARNING: matrix D is non-zero. Algebraic loops may occur.')
        end
    end
end

object.controller = controller;

xnames = controller.getStateNames();
unames = controller.getInputNames();
pnames = controller.getParameterNames();
dnames = controller.getUnmeasurableInputNames();

disp('The controller has been associated to the system.')
disp('Check the names below to ensure correct correspondence:')
for i = 1:object.nx
    disp([object.xnames{i},' <---> ',xnames{i}])    
end
for i = 1:object.nu
    disp([object.unames{i},' <---> ',unames{i}])    
end
for i = 1:object.np
    disp([object.pnames{i},' <---> ',pnames{i}])    
end
for i = 1:object.nd
    disp([object.dnames{i},' <---> ',dnames{i}])    
end
disp(' ')