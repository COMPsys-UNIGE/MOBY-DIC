classdef embeddedSystem
    % embeddedSystem   embedded system used to observe and control a LTI or
    %                  PWA system
    %
    % OBJ = embeddedSystem()
    % Builds an empty embeddedSystem object OBJ.
    %
    % OBJ = embeddedSystem(dynSys)
    % Builds an embeddedSystem object from a dynSys object, retrieving the 
    % observer and the controller associated to it.
    %
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
    
    
    properties (Access = protected)
        nx = 0; % Number of states
        nu = 0; % Number of inputs
        np = 0; % Number of parameters
        nd = 0; % Number of unmeasurable inputs
        ny = 0; % Number of outputs
        xnames = [];    % Names of the states
        unames = [];    % Names of the inputs
        pnames = [];    % Names of the parameters
        dnames = [];    % Names of the unmeasurable inputs
        ynames = [];    % Names of the outputs
        controller = [];    % Controller for the dynamical system
        observer = [];    % Observer for the dynamical system
%         constraints = [];   % Linear inequality constraints
        
        range = [];
        
        dynSys = [];
        
    end
    
    
    % Methods
    
    methods
        
        function es = embeddedSystem(varargin)
            if nargin == 1 || nargin == 2
                if ~isa(varargin{1},'dynSys')
                    error('Input must be a dynSys object');
                end
                sis = varargin{1};
                es.dynSys = sis;
                es.nx = sis.getNumberOfStates;
                es.nu = sis.getNumberOfInputs;
                es.ny = sis.getNumberOfOutputs;
                es.np = sis.getNumberOfParameters;
                es.nd = sis.getNumberOfUnmeasurableInputs;
                
                if es.nx <= 0
                    error('Number of states must be > 0')
                end
                if es.nu < 0
                    error('Number of inputs must be >= 0')
                end
                if es.ny <= 0
                    error('Number of outputs must be > 0')
                end
                if es.np < 0
                    error('Number of parameters must be >= 0')
                end
                if es.nd < 0
                    error('Number of unmeasurable inputs must be >= 0')
                end
                
                
                es.xnames = sis.getStateNames;
                
                es.unames = sis.getInputNames;
                
                es.ynames = sis.getOutputNames;
                
                es.pnames = sis.getParameterNames;
                
                es.dnames = sis.getUnmeasurableInputNames;
                
                es.controller = sis.getController;
                
                es.observer = sis.getObserver;
                
%                 es.constraints = [];
                
                if nargin == 2
                    r = varargin{2};
                else
                    r = [];
                end
                
                % check over range input
                if es.hasController && ~es.hasObserver
                    if ~isempty(r)
                        warning('Embedded system use controller parameter instead of parameters in range structure');
                    end
                    
%                     es.constraints = es.controller.getConstraints;
                    
                    es.range = es.getRangeFromController();
                    
                elseif es.hasObserver && ~es.hasController
                    if ~isstruct(r)
                        error('range must be a struct');
                    else
                        if ~isfield(r,'xmin')||~isfield(r,'xmax')||...
                                ~isfield(r,'umin')||~isfield(r,'umax')||...
                                ~isfield(r,'pmin')||~isfield(r,'pmax')||...
                                ~isfield(r,'dmin')||~isfield(r,'dmax')||...
                                ~isfield(r,'ymin')||~isfield(r,'ymax')
                            error(['range options must be a struct with fields xmin,xmax,umin,umax,'...
                                'pmin,pmax,dmin,dmax,ymin,ymax']);
                        else
                            if numel(r.xmin) == 1
                                es.range.xmin = repmat(r.xmin,es.nx);
                                es.range.xmax = repmat(r.xmax,es.nx);
                            else
                                es.range.xmin = r.xmin;
                                es.range.xmax = r.xmax;
                            end
                            if numel(r.umin) == 1
                                es.range.umin = repmat(r.umin,es.nu);
                                es.range.umax = repmat(r.umax,es.nu);
                            else
                                es.range.umin = r.umin;
                                es.range.umax = r.umax;
                            end
                            if numel(r.dmin) == 1
                                es.range.dmin = repmat(r.dmin,es.nd);
                                es.range.dmax = repmat(r.dmax,es.nd);
                            else
                                es.range.dmin = r.dmin;
                                es.range.dmax = r.dmax;
                            end
                            if numel(r.pmin) == 1
                                es.range.pmin = repmat(r.pmin,es.np);
                                es.range.pmax = repmat(r.pmax,es.np);
                            else
                                es.range.pmin = r.pmin;
                                es.range.pmax = r.pmax;
                            end
                            if numel(r.ymin) == 1
                                es.range.ymin = repmat(r.ymin,es.ny);
                                es.range.ymax = repmat(r.ymax,es.ny);
                            else
                                es.range.ymin = r.ymin;
                                es.range.ymax = r.ymax;
                            end
                            
                        end
                    end
                elseif es.hasObserver && es.hasController
%                     es.constraints = es.controller.getConstraints;

                    if isa(es.getController, 'implicitMPCctrl')
                        if ~isstruct(r)
                            error('range must be a struct');
                        else
                            if ~isfield(r,'xmin')||~isfield(r,'xmax')||...
                                    ~isfield(r,'umin')||~isfield(r,'umax')||...
                                    ~isfield(r,'pmin')||~isfield(r,'pmax')||...
                                    ~isfield(r,'dmin')||~isfield(r,'dmax')||...
                                    ~isfield(r,'ymin')||~isfield(r,'ymax')
                                error(['range options must be a struct with fields xmin,xmax,umin,umax,'...
                                    'pmin,pmax,dmin,dmax,ymin,ymax']);
                            else
                                if numel(r.xmin) == 1
                                    es.range.xmin = repmat(r.xmin,es.nx);
                                    es.range.xmax = repmat(r.xmax,es.nx);
                                else
                                    es.range.xmin = r.xmin;
                                    es.range.xmax = r.xmax;
                                end
                                if numel(r.umin) == 1
                                    es.range.umin = repmat(r.umin,es.nu);
                                    es.range.umax = repmat(r.umax,es.nu);
                                else
                                    es.range.umin = r.umin;
                                    es.range.umax = r.umax;
                                end
                                if numel(r.dmin) == 1
                                    es.range.dmin = repmat(r.dmin,es.nd);
                                    es.range.dmax = repmat(r.dmax,es.nd);
                                else
                                    es.range.dmin = r.dmin;
                                    es.range.dmax = r.dmax;
                                end
                                if numel(r.pmin) == 1
                                    es.range.pmin = repmat(r.pmin,es.np);
                                    es.range.pmax = repmat(r.pmax,es.np);
                                else
                                    es.range.pmin = r.pmin;
                                    es.range.pmax = r.pmax;
                                end
                                if numel(r.ymin) == 1
                                    es.range.ymin = repmat(r.ymin,es.ny);
                                    es.range.ymax = repmat(r.ymax,es.ny);
                                else
                                    es.range.ymin = r.ymin;
                                    es.range.ymax = r.ymax;
                                end
                            end
                        end
                        es.range = r;
                    else
                        es.range = es.getRangeFromController();
                    end
                    
                    if ~isempty(r)
                        if ~isstruct(r) || ~isfield(r,'ymin') || ~isfield(r,'ymax')
                            error('range must be a struct with fields ymin and ymax')
                        else
                            if numel(r.ymin) == es.getNumberOfOutputs && numel(r.ymax) == es.getNumberOfOutputs
                                es.range.ymin = r.ymin;
                                es.range.ymax = r.ymax;
                            elseif numel(r.ymin) == 1 && numel(r.ymax) == 1
                                es.range.ymin = repmat(r.ymin,es.getNumberOfOutputs,1);
                                es.range.ymax = repmat(r.ymax,es.getNumberOfOutputs,1);
                            else
                                error('fields ymin and ymax of struct range must have lenght ny or 1');
                            end
                        end
                    else
                        error('If observer is provided output variation range must be provided.');
                    end
                else
                    error('In order to create embedded system dynSys must have a controller or an observer');
                end
            end
        end
        
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
        
        generateC(object,varargin);
        
        varargout = generateVHDL(object,circuit_parameters,options);
        
        disp(object);
        
    end
    
end
