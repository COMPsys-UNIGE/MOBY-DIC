classdef observer
    
    % observer   Generic state observer for a dynSys object
    %
    % This object represents a generic state observer for a dynamical 
    % system (dynSys object).
    %
    % This is an abstract class which cannot be istantiated.
    %
    % observer methods:
    % getInputNames - gets the mnemonical name of the input variables
    % getNumberOfUnmeasurableInputs - gets the number of unmeasurable inputs
    % getNumberOfInputs - gets the number of input variables (control actions)
    % getNumberOfOutputs - gets the number of system outputs
    % getNumberOfParameters - gets the number of parameters
    % getNumberOfStates - gets the number of state variables
    % getOutputNames - gets the mnemonical name of the output variables
    % getParameterNames - gets the mnemonical name of the parameters
    % getSamplingTime - gets the sampling time of the system
    % getStateNames - gets the mnemonical name of the state variables
    % getUnmeasurableInputNames - gets the mnemonical name of the unmeasurable inputs
    % setInputNames - gives a mnemonical name to the input variables (control actions)
    % getNumberOfOutputs - gives a mnemonical name to the system outputs
    % setParameterNames - gives a mnemonical name to the parameters
    % setStateNames - gives a mnemonical name to the state variables
    % setUnmeasurableInputNames - gives a mnemonical name to the unmeasurable inputs
    %
    % See also kalmanFilter, kalmanPredictor.
    
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
      
    % Properties
    
    % TO DO
    % Mettere metodi predict e update e computeOutput come virtuali
    
    properties (Access = protected)
        nx = 0; % Number of states
        nu = 0; % Number of inputs
        np = 0; % Number of parameters
        nd = 0; % Number of unmeasurable inputs
        ny = 0; % Number of outputs
        Ts = 0;         % Sampling time
        xnames = [];    % Names of the states
        unames = [];    % Names of the inputs
        pnames = [];    % Names of the parameters
        ynames = [];    % Names of system outputs
        dnames = [];    % Names of the unmeasurable inputs
    end
    
    methods (Abstract)
                       
        % Other methods
                       
        disp(object);
        
        xpred = predict(object,xcur,ucur,pcur,ycur);
        
        [xcur, ycur] = computeOutput(object,x,u,p,y);
        
    end
    
    methods
        
        function obs = observer()
            % abstract class empty cconstructor
        end
                
        % Get methods
        
        nx = getNumberOfStates(object);
        
        nu = getNumberOfInputs(object);
        
        ny = getNumberOfOutputs(object);
        
        np = getNumberOfParameters(object);
        
        nd = getNumberOfUnmeasurableInputs(object);
        
        Ts = getSamplingTime(object);
        
        names = getStateNames(object);
        
        names = getInputNames(object);
        
        names = getParameterNames(object);
        
        names = getUnmeasurableInputNames(object);
                
        % Set methods
        
        object = setStateNames(object,stateNames);
        
        object = setInputNames(object,inputNames);
        
        object = setOutputNames(object,outputNames);
        
        object = setParameterNames(object,parameterNames);
        
        object = setUnmeasurableInputNames(object,unmeasurableInputNames);

        
    end
end