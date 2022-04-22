classdef constraints
    % constraints   Object describing linear inequality constraints on states,
    %               parameters, outputs, unmeasurable and measurable inputs 
    %               of a dynamical system
    %
    % This obect represents linear inequality constraints involving states
    % (x), parameters (p), unmeasurable (d) and measurable (u) inputs and 
    % outputs (y) of a dynamical system as well as the reference state or 
    % output (ref) in case of tracking control. The constraints are in the 
    % form:
    %    _   _
    %   |  x  |
    %   |  u  | 
    % H |  p  | <= K
    %   |  d  |
    %   |  y  |
    %   |_ref_|
    %
    % Constraints can be set for different time instants, up to a time
    % horizon N. Either hard or soft constraints can be defined.
    %
    % OBJ = constraints()
    % Builds an empty constraints object OBJ.
    %
    % OBJ = constraints(NX,NU,NY,N)
    % Builds a constraints object OBJ by defining the number of state
    % variables (NX), the number of input variables (NU), the number of
    % outputs (NY) and the time horizon (N). The number of parameters, 
    % unmeasurable inputs and reference states is assumed to be 0.
    %
    % OBJ = constraints(NX,NU,NY,NREF,N)
    % Builds a constraints object OBJ by defining the number of state
    % variables (NX), the number of input variables (NU), the number of
    % outputs (NY), the number of reference states or outputs (NREF) and 
    % the time horizon (N). The number of parameters and unmeasurable 
    % outputs is assumed to be 0.
    %
    % OBJ = constraints(NX,NU,NY,NP,ND,N)
    %
    % Builds a constraints object OBJ by defining the number of state
    % variables (NX), the number of input variables (NU), the number of
    % outputs (NY), the number of parameters (NP), the number of 
    % unmeasurable inputs (ND) and the time horizon (N). The number of 
    % reference states is assumed to be 0.
    %
    % OBJ = constraints(NX,NU,NY,NP,ND,NREF,N)
    %
    % Builds a constraints object OBJ by defining the number of state
    % variables (NX), the number of input variables (NU), the number of
    % outputs (NY), the number of parameters (NP), the number of unmeasurable 
    % inputs (ND), the number of reference states or outputs (NREF) and the time 
    % horizon (N). The reference states are meaningful only for tracking control.
    %
    % constraints methods:
    % disp - displays some information about the constraints object.
    % getAllConstraints - gets the H and K matrices defining all constraints.
    % getInputConstraints - gets the H and K matrices defining input constraints.
    % getNumberOfInputs - gets the number of input variables.
    % getNumberOfOutputs - gets the number of output variables.
    % getNumberOfParameters - gets the number of parameters.
    % getNumberOfReferences - gets the number of reference states or outputs.
    % getNumberOfStates - gets the number of state variables.
    % getNumberOfUnmeasurableInputs - gets the number of unmeasurable inputs.
    % getOutputConstraints - gets the H and K matrices defining output constraints.
    % getParameterConstraints - gets the H and K matrices defining parameter constraints.
    % getReferenceConstraints - gets the H and K matrices defining reference constraints.
    % getStateConstraints - gets the H and K matrices defining state constraints.
    % getTimeHorizon - gets the time horizon.
    % getUnmeasurableInputConstraints - gets the H and K matrices defining unmeasurable input constraints.
    % setConstraints - sets inequality constraints.
    % setSoftConstraints - sets soft inequality constraints.
    %
    % See also MPCctrl.
    
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
    
    properties (Access = private)
        nx = 0;
        nu = 0;
        np = 0;
        nd = 0;
        ny = 0;
        nref = 0;
        N = 0;
        H = [];
        K = [];
        soft = [];
    end
    
    % Methods
    
    methods
        
        % Constructor
        function constr = constraints(varargin)
            
            % Empty object
            if nargin == 0
                constr.N = 0;
                constr.H = [];
                constr.K = [];
                constr.soft = [];
                constr.nx = 0;
                constr.nu = 0;
                constr.np = 0;
                constr.nd = 0;
                constr.ny = 0;
                constr.nref = 0;
                
            elseif nargin == 4
                constr.nx = varargin{1};
                constr.nu = varargin{2};
                constr.ny = varargin{3};
                constr.np = 0;
                constr.nd = 0;
                constr.nref = 0;
                constr.N = varargin{4};
                constr.H = cell(constr.N+1,1);
                constr.K = cell(constr.N+1,1);
                constr.soft = cell(constr.N+1,1);
                
                if constr.nx <= 0
                    error('Number of states must be at least 1.')
                end
                if constr.N <= 0
                    error('Prediction horizon must be at least 1.')
                end
                if constr.nu < 0
                    error('nu must be positive.')
                end
                if constr.ny < 0
                    error('ny must be positive.')
                end
                
            elseif nargin == 5
                constr.nx = varargin{1};
                constr.nu = varargin{2};
                constr.ny = varargin{3};
                constr.np = 0;
                constr.nd = 0;
                constr.nref = varargin{4};
                constr.N = varargin{5};
                constr.H = cell(constr.N+1,1);
                constr.K = cell(constr.N+1,1);
                constr.soft = cell(constr.N+1,1);
                
                if constr.nx <= 0
                    error('Number of states must be at least 1.')
                end
                if constr.N <= 0
                    error('Prediction horizon must be at least 1.')
                end
                if constr.nu < 0
                    error('nu must be positive.')
                end                
                if constr.nref < 0
                    error('nref must be positive.')
                end
                if constr.ny < 0
                    error('ny must be positive.')
                end
                
            elseif nargin == 6
                constr.nx = varargin{1};
                constr.nu = varargin{2};
                constr.ny = varargin{3};
                constr.np = varargin{4};
                constr.nd = varargin{5};
                constr.nref = 0;
                constr.N = varargin{6};
                constr.H = cell(constr.N+1,1);
                constr.K = cell(constr.N+1,1);
                constr.soft = cell(constr.N+1,1);
                
                if constr.nx <= 0
                    error('Number of states must be at least 1.')
                end
                if constr.N <= 0
                    error('Prediction horizon must be at least 1.')
                end
                if constr.nu < 0
                    error('nu must be positive.')
                end
                if constr.np < 0
                    error('np must be positive.')
                end
                if constr.nd < 0
                    error('nd must be positive.')
                end
                if constr.ny < 0
                    error('ny must be positive.')
                end
                
            elseif nargin == 7
                constr.nx = varargin{1};
                constr.nu = varargin{2};
                constr.ny = varargin{3};
                constr.np = varargin{4};
                constr.nd = varargin{5};
                constr.nref = varargin{6};
                constr.N = varargin{7};
                constr.H = cell(constr.N+1,1);
                constr.K = cell(constr.N+1,1);
                constr.soft = cell(constr.N+1,1);
                
                if constr.nx <= 0
                    error('Number of states must be at least 1.')
                end
                if constr.N <= 0
                    error('Prediction horizon must be at least 1.')
                end
                if constr.nu < 0
                    error('nu must be positive.')
                end
                if constr.np < 0
                    error('np must be positive.')
                end
                if constr.nd < 0
                    error('nd must be positive.')
                end
                if constr.nref < 0
                    error('nref must be positive.')
                end
                if constr.ny < 0
                    error('ny must be positive.')
                end
                
            else
                error('Wrong input arguments for constraints object constructor');
                
            end
        end
        
        % Set methods
        
        object = setConstraints(object,varargin);
        object = setSoftConstraints(object,varargin);
        
        % Get methods
        
        nx = getNumberOfStates(object);
        
        ny = getNumberOfOutputs(object);
        
        nref = getNumberOfReferences(object);
        
        nu = getNumberOfInputs(object);
        
        nu = getNumberOfParameters(object);
        
        nd = getNumberOfUnmeasurableInputs(object);
        
        N = getTimeHorizon(object);
        
        [H, K] = getAllConstraints(object,varargin);
        
        [H, K] = getStateConstraints(object,varargin);
        
        [H, K] = getOutputConstraints(object,varargin);
        
        [H, K] = getReferenceConstraints(object,varargin);
        
        [H, K] = getInputConstraints(object,varargin);
        
        [H, K] = getParameterConstraints(object,varargin);
        
        [H, K] = getUnmeasurableInputConstraints(object,varargin);
        
        % Other methods
        
        disp(object);
        
    end
end


