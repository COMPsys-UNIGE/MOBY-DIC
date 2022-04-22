classdef controller
      
    % controller   Generic controller for a dynSys object
    %
    % This object represents a generic controller for the regulation or
    % tracking of a dynamical system (dynSys object).
    %
    % This is an abstract class which cannot be istantiated.
    %
    % controller methods:
    % getControlVariable - indicates if the controller regulates (or tracks)
    %                      system states or outputs
    % getDomain - gets the controller domain
    % getFunction - gets the control function
    % getInputNames - gets the mnemonical name of the input variables
    % getNumberOfDimensions - gets the number of dimensions of the control function 
    % getNumberOfInputs - gets the number of input variables (control actions)
    % getNumberOfOutputs - gets the number of output variables
    % getNumberOfParameters - gets the number of parameters
    % getNumberOfStates - gets the number of state variables
    % getNumberOfUnmeasurableInputs - gets the number of unmeasurable inputs
    % getOutputNames - gets the mnemonical name of the outputs
    % getParameterNames - gets the mnemonical name of the parameters
    % getReference - gets the reference state or output
    % getSamplingTime - gets the sampling time of the system
    % getStateNames - gets the mnemonical name of the state variables
    % getTrackingVariable - gets the index of the tracking variable
    % getUnmeasurableInputNames - gets the mnemonical name of the unmeasurable inputs
    % isTracking - Returns true if the control function is for tracking, false if it is for regulation    
    % setInputNames - gives a mnemonical name to the input variables (control actions)
    % setOutputNames - gives a mnemonical name to the output variables 
    % setParameterNames - gives a mnemonical name to the parameters
    % setStateNames - gives a mnemonical name to the state variables
    % setUnmeasurableInputNames - gives a mnemonical name to the unmeasurable inputs
    %
    % See also MPCctrl
    
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
        nfun = 0;   % Number of control functions (i.e., one for each time step)
        tracking = false;   % Tracking or regulation controller
        trackvar = [];   % Index of variable to track
        ctrlvar = '';   % Indicates if the controller regulates (or tracks) states or outputs
        ref = [];      % Reference state or output (only for regulation)
        Ts = 0;         % Sampling time
        xnames = [];    % Names of the states
        unames = [];    % Names of the inputs
        pnames = [];    % Names of the parameters
        dnames = [];    % Names of the unmeasurable inputs
        ynames = [];    % Names of the outputs
        fun = [];       % Control function
    end
    
    methods (Abstract)
                       
        % Other methods
        
        varargout = eval(object,x,varargin);
        
        generateC(object,varargin);
                
        disp(object);
        
    end
    
    methods
        
        function ctrl = controller()
            % abstract class empty constructor
        end
                
        % Get methods
        
        ndim = getNumberOfDimensions(object);
        
        nx = getNumberOfStates(object);
        
        nu = getNumberOfInputs(object);
        
        ny = getNumberOfOutputs(object);
        
        np = getNumberOfParameters(object);
        
        nd = getNumberOfUnmeasurableInputs(object);
        
        ref = getReference(object);
        
        Ts = getSamplingTime(object);
        
        names = getStateNames(object);
        
        names = getInputNames(object);
        
        names = getOutputNames(object);
        
        names = getParameterNames(object);
        
        names = getUnmeasurableInputNames(object);
        
        trackvar = getTrackingVariable(object);
        
        ctrlvar = getControlVariable(object);
        
        [Hd, Kd] = getDomain(object);
        
        func = getFunction(object);
                
        % Is methods
        
        answ = isTracking(object);
        
        % Set methods
        
        object = setStateNames(object,stateNames);
        
        object = setInputNames(object,inputNames);
        
        object = setOutputNames(object,outputNames);
        
        object = setParameterNames(object,parameterNames);
        
        object = setUnmeasurableInputNames(object,unmeasurableInputNames);

        
    end
end