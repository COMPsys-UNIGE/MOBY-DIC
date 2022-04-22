function [H, K] = getAllConstraints(object,varargin)
% getAllConstraints   Gets the H and K matrices defining all linear inequality constraints
%
% [H, K] = getAllConstraints(OBJ)
% Gets all constraints. H and K are cell arrays, whose elements correspond
% to different time instants (up to time horizon N). Each element of the
% cell arays defines a system of inequalities in the form:
%    _   _
%   |  x  |
%   |  u  |
% H |  p  | <= K
%   |  d  |
%   |  y  |
%   |_ref_|
%
% being x, u, p, d, y and ref the system states, inputs, parameters, 
% unmeasurable inputs, outputs and references, respectively. 
% OBJ is the constraints object.
%
% [H, K] = getAllConstraints(OBJ,T)
% Gets all constraints at time instant T. T must be a scalar between 0 and
% the time horizon N. H and K are matrices defining all constraints in the
% form:
%    _   _
%   |  x  |
%   |  u  |
% H |  p  | <= K
%   |  d  |
%   |  y  |
%   |_ref_|
%
% [H, K] = getAllConstraints(OBJ,'soft')
% [H, K] = getAllConstraints(OBJ,'hard')
% Gets only soft or hard constraints.
%
% [H, K] = getAllConstraints(OBJ,'soft',T)
% [H, K] = getAllConstraints(OBJ,'hard',T)
% Gets only soft or hard constraints at time instant T.
%
% See also: constraints/getStateConstraints,
% constraints/getInputConstraints, constraints/getParameterConstraints, 
% constraints/getUnmeasurableInputConstraints

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



% Extract all matrices
if nargin == 1
    
    H = object.H;
    K = object.K;
    
    
elseif nargin == 2
    
    % Extract matrices at a given time instant
    if isnumeric(varargin{1})
        
        T = varargin{1};
        
        if numel(T) > 1
            error('T must be a scalar')
        end
        
        if T < 0 || T > object.N
            error(['time must be comprised between 0 and ',num2str(object.N),'!'])
        end
        
        H = object.H{T+1};
        K = object.K{T+1};
        
    elseif ischar(varargin{1})
        
        label = varargin{1};
        
        % Extracting only soft constraints
        if strcmpi(label,'soft')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx = object.soft{i} == 1;
                
                H{i} = object.H{i}(idx,:);
                K{i} = object.K{i}(idx);
                
            end
            
            % Extracting only hard constraints
        elseif strcmpi(label,'hard')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx = object.soft{i} == 0;
                
                H{i} = object.H{i}(idx,:);
                K{i} = object.K{i}(idx);
                
            end
            
        else
            
            error('Label must be either ''soft'' or ''hard''');
            
        end
        
    else
        error('Wrong input arguments');
    end

elseif nargin == 3
    
    label = varargin{1};
    T = varargin{2};
    
    if numel(T) > 1
        error('T must be a scalar')
    end
    
    if T < 0 || T > object.N
        error(['time must be comprised between 0 and ',num2str(object.N),'!'])
    end
    
    % Extracting only soft constraints at time T
    if strcmpi(label,'soft')
        
        idx = object.soft{T+1} == 1;
    
        H = object.H{T+1}(idx,:);
        K = object.K{T+1}(idx);
        
        % Extracting only hard constraints at time T
    elseif strcmpi(label,'hard')
        
        idx = object.soft{T+1} == 0;
    
        H = object.H{T+1}(idx,:);
        K = object.K{T+1}(idx);
        
    else
        
        error('Label must be either ''soft'' or ''hard''');
        
    end
    
else
    
    error('Wrong input number');
end
