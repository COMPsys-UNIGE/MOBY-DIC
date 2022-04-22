classdef discontinuousPwasFunction < MOBYFunction    
    % discontinuous pwasFunction DiscontinuousPiece-Wise Affine 
    %                            Simplicial function
    % This object represents a discontinuous piecewise affine function 
    % defined over a uniform or non-uniform simplicial domain partition.
      
    % OBJ = discontinuousPwasFunction()
    % Builds an empty discontinuous_pwasFunction object OBJ.
    %
    % OBJ = discontinuousPwasFunction(FUN,H,K)
    % Builds a discontinuous pwasFunction object OBJ by specifying the
    % array of pwas function FUN and the arrays H/K that contain the  
    % specified constrains for each domain.
 
    % discontinuousPwasFunction methods:
    %   disp - displays some information about the discontinuous pwasFunction object.
    %   eval - evaluates the discontinuous PWAS function.
    %   plot - plots the discontinuous PWAS function.
    %   findDynamics - find the index of the dynamic.
    %
    % The discontinuous pwasFunction object is derived from MOBYFunction 
    % and inherits all its methods.
    %
    % See also MOBYFunction, pwasFunction.
    
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
    %---------------------------------------------------------------------
        
    properties (Access = private)
        
        nDyn = 0;
        functions = [];
        
    end
    
    methods
        
        function  object = discontinuousPwasFunction(varargin)
            
            if nargin == 0
                
                object.nDyn = 0;
                object.functions = [];
                
            elseif nargin == 3   % object=discontinuous_pwasFunction (fun,H,K)
                
                % Extract inputs
                fun = varargin{1};
                Hmat = varargin {2};
                Kmat = varargin {3};
                
                % Check inputs
                nfun = numel (fun);
                object.nDyn = nfun;
                nH = numel (Hmat);
                nK = numel (Kmat);
                
                if nfun < 2
                    error('A discontinuous pwas funcion needs at least 2 dynamics');
                end
                if nH ~= nK
                    error('The number of K matrices must be the same of the number of the K matrices');
                elseif nH ~= object.nDyn
                    error('The number of dynamics must be the same of the number of the K/H matrices');
                end
                
                for i = 1:object.nDyn
                    if ~isa(fun{i},'pwasFunction')
                        error('First argument must be pwas function object');
                    end
                    object.functions(i).H = Hmat{i};
                    object.functions(i).K = Kmat{i};
                    object.functions(i).pwasFunction = fun{i};
                end
                
                object.nx = fun{1}.getDomainDimensions;
                object.ny = fun{1}.getCodomainDimensions;
                object.xnames = fun{1}.getInputNames;
                object.ynames = fun{1}.getOutputNames;
                [Hd Kd] = fun{1}.getDomain;
                object.domain.Hd = Hd;
                object.domain.Kd = Kd;
                
            end
            
        end
        
        plot(object,varargin);
        
        y = eval(object,x,options);
        
        disp(object);
        
        dyn = findDynamics(object,x);
        
        [Hd, Kd] = getDomain(object);
        
        uniform = isUniform(object);
        
        np = getNumberOfPartitions(object);
        
        P = getPartition(object);
        
        V = getVertices(object,varargin)

    end
    
end

