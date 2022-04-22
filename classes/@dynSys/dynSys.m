classdef dynSys
    
    % dynSys   Generic dynamical system
    %
    % This object represents a generic dynamical system.
    % The system is characterized by a set of states (x), inputs (u),
    % measurable parameters (p), unmeasurable inputs (d) and outputs (y).
    %
    % The continuous time formulation is:
    %       dx/dt = f(x,u,p,d)
    %           y = h(x,u,p,d)
    %
    % while the discrete time formulation is:
    %      x(k+1) = f(x(k),u(k),p(k),d(k))
    %        y(k) = h(x(k),u(k),p(k),d(k))
    %
    % This is an abstract class which cannot be istantiated.
    %
    % dynSys methods:
    %   getController - gets the controller associated to the system    
    %   getInputNames - gets the mnemonical name of the input    
    %   getNumberOfInputs - gets the number of inputs
    %   getNumberOfOutputs - gets the number of outputs
    %   getNumberOfParameters - gets the number of parameters
    %   getNumberOfStates - gets the number of states
    %   getNumberOfUnmeasurableInputs - gets the number of unmeasurable inputs
    %   getObserver - gets the observer associated to the system
    %   getOutputNames - gets the mnemonical name of the outputs
    %   getParameterNames - gets the mnemonical name of the parameters
    %   getSamplingTime - gets the sampling time of the system
    %   getStateNames - gets the mnemonical name of the states
    %   getUnmeasurableInputNames - gets the mnemonical name of the unmeasurable inputs
    %   hasController - says if a controller is associated to this system
    %   hasObserver - says if an observer is associated to this system
    %   isContinuousTime - indicates if the system is continuous-time
    %   isDiscreteTime - indicates if the system is discrete-time
    %   setController - associates a controller to the system
    %   setInputNames - gives a mnemonical name to the inputs
    %   setObserver - associates an observer to this system
    %   setOutputNames - gives a mnemonical name to the outputs
    %   setParameterNames - gives a mnemonical name to the parameters
    %   setStateNames - gives a mnemonical name to the states
    %   setUnmeasurableInputNames - gives a mnemonical name to the unmeasurable inputs
    %
    % See also ltiSys, pwaSys
    
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
    

    % Properties
    
    properties (Access = protected)
        nx = 0; % Number of states
        nu = 0; % Number of inputs
        np = 0; % Number of parameters
        nd = 0; % Number of unmeasurable inputs
        ny = 0; % Number of outputs
        mode = ''; % Discrete- or continuous-time
        Ts = 0;         % Sampling time
        xnames = [];    % Names of the states
        unames = [];    % Names of the inputs
        pnames = [];    % Names of the parameters
        dnames = [];    % Names of the unmeasurable inputs
        ynames = [];    % Names of the outputs
        controller = [];    % Controller for the dynamical system
        observer = [];    % Observer for the dynamical system
    end
    
    methods (Abstract)
                
        % Other methods
        
        y = computeOutput(object,x,u,p,d)
        
        xnext = evaluate(object,x,u,p,d)
        
        [X, U, Y, T] = sim(object,mode,x0,p,d,N,ctrl);
        
        [X, U, Y, T] = simplot(object,mode,x0,p,d,N,ctrl);
        
        dynSys = discretize(object,Ts);
        
        disp(object);
        
    end
    
    methods
        
        function sys = dynSys()
            % abstract class empty cconstructor
        end
        
        % Is methods
        
        ctime = isContinuousTime(object);
        dtime = isDiscreteTime(object);
        
        % Has methods
        
        answ = hasController(object);
        
        answ = hasObserver(object);
        
        % Get methods
        
        nx = getNumberOfStates(object);
        
        nu = getNumberOfInputs(object);
        
        ny = getNumberOfOutputs(object);
        
        np = getNumberOfParameters(object);
        
        nd = getNumberOfUnmeasurableInputs(object);
        
        Ts = getSamplingTime(object);
        
        observer = getObserver(object);
        
        controller = getController(object);
               
        names = getStateNames(object);
        
        names = getInputNames(object);
        
        names = getOutputNames(object);
        
        names = getParameterNames(object);
        
        names = getUnmeasurableInputNames(object);
        
        % Set methods
        
        object = setStateNames(object,stateNames);
        
        object = setInputNames(object,inputNames);
        
        object = setOutputNames(object,outputNames);
        
        object = setParameterNames(object,parameterNames);
        
        object = setUnmeasurableInputNames(object,unmeasurableInputNames);
        
        object = setObserver(object,observer);
        
        object = setController(object,controller);
                
    end
end