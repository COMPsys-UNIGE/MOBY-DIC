classdef MOBYFunction
    
    % MOBYFunction   Generic function object of MOBY-DIC Toolbox
    %
    % This object represents a generic vector function f defined over a 
    % polytopic domain.
    % This is an abstract class which cannot be istantiated.
    %
    % MOBYFunction methods:
    %   approximate - approximates the function with a PWAS function
    %   getCodomainDimensions - gets the number of codomain dimensions    
    %   getDomain - gets the edges of the function domain
    %   getDomainDimensions - gets the number of domain dimensions
    %   getInputNames - gets the names of the function inputs
    %   getOutputNames - gets the names of the function outputs
    %   setInputNames - sets the names of the function inputs
    %   setOutputNames - sets the names of the function outputs
    %
    % See also pwagFunction, pwasFunction, nonlinFunction
    
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
        nx = 0; % Number of domain dimensions
        ny = 0; % Number of codomain dimensions
        domain = struct('Hd',[],'Kd',[]);   % Struct defining function domain
        xnames = [];    % Input names
        ynames = [];    % Output names
    end
    
    
    % Methods
    methods (Abstract)
        
        plot(object,varargin);
        
        y = eval(object,x,options);
        
        disp(object);
        
    end
    
    methods
        
        fappr = approximate(object,options);
        
        [Hd, Kd] = getDomain(object);
        
        nx = getDomainDimensions(object);
        
        ny = getCodomainDimensions(object);
        
        xnames = getInputNames(object);
        
        ynames = getOutputNames(object);
        
        object = setInputNames(object,xnames);
        
        object = setOutputNames(object,ynames);
        
    end
end


