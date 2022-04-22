function object = setObserver(object,observer)
% setObserver   Associates an observer to the system
%
% OBJ = setObserver(OBJ,OBS)
% OBS must be an observer object. OBJ is the dynSys object.
%
% The observer must have the same number of inputs, states, parameters,
% outputs and unmeasurable inputs. If the system is discrete-time, the
% sampling time of the observer must be equal or multiple of the system sampling
% time. If a controller is already associated to the system, the sampling
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


% Check compatibility

if ~isa(observer,'observer')
    error('OBS must be a observer object')
end

if observer.getNumberOfStates() ~= object.nx
    error('Observer and system have different number of states');
end
if observer.getNumberOfInputs() ~= object.nu
    error('Observer and system have different number of inputs');
end
if observer.getNumberOfOutputs() ~= object.ny
    error('Observer and system have different number of outputs');
end
if observer.getNumberOfParameters() ~= object.np
    error('Observer and system have different number of parameters');
end
if observer.getNumberOfUnmeasurableInputs() ~= object.nd
    error('Observer and system have different number of unmeasurable inputs');
end

% Retrieve rounding tolerance
pars = getMOBYDICpars();
tol = pars.roundtol;

% Check sampling time
if object.isDiscreteTime()
    if observer.getSamplingTime() < object.Ts
        error('Sampling time of the observer must not be lower than system sampling time');
    end
    
    % Ratio between sampling times
    Tsratio = observer.getSamplingTime()/object.Ts;
    
    % Round sampling time ratio, to prevent numerical issues
    Tsratio = round(Tsratio/tol)*tol;
    
    if floor(Tsratio) ~= Tsratio
        error('The sampling time of the observer must be a multiple of the system sampling time');
    end
end

% Check sampling time
if object.hasController()
    ctrl = object.controller;
    if ctrl.getSamplingTime() < observer.getSamplingTime()
        error('Sampling time of the controller must not be lower than observer sampling time');
    end
    
    % Ratio between sampling times
    Tsratio = ctrl.getSamplingTime()/observer.getSamplingTime();
    
    % Round sampling time ratio, to prevent numerical issues
    Tsratio = round(Tsratio/tol)*tol;
    
    if floor(Tsratio) ~= Tsratio
        error('The sampling time of the controller must be a multiple of the observer sampling time');
    end
    
    if isa(observer,'kalmanFilter') && isa(object,'ltiSys')
        D = object.getMatrices('D');
        if any(D ~= 0)
            disp('WARNING: matrix D is non-zero. Algebraic loops may occur.')
        end
    end
end

object.observer = observer;

xnames = observer.getStateNames();
unames = observer.getInputNames();
pnames = observer.getParameterNames();
dnames = observer.getUnmeasurableInputNames();

disp('The observer has been associated to the system.')
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

end