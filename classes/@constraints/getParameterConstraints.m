function [H, K] = getParameterConstraints(object,varargin)
% getParameterConstraints   Gets the H and K matrices defining parameter constraints
%
% [H, K] = getParameterConstraints(OBJ)
% Gets all constraints involving system parameters. H and K are cell arrays, 
% whose elements corresponds to different time instants (up to time horizon N). 
% Each elements of the cell arays defines a system of inequalities in the form:
%
% H p  <= K
%
% being p the system parameters. OBJ is the constraints object.
%
% [H, K] = getParameterConstraints(OBJ,T)
% Gets all constraints involving system parameters at time instant T. 
% T must be a scalar between 0 and the time horizon N. 
% H and K are matrices defining all parameter constraints in the form:
%
% H p  <= K
%
% [H, K] = getParameterConstraints(OBJ,'soft')
% [H, K] = getParameterConstraints(OBJ,'hard')
% Gets only soft or hard parameter constraints.
%
% [H, K] = getParameterConstraints(OBJ,'soft',T)
% [H, K] = getParameterConstraints(OBJ,'hard',T)
% Gets only soft or hard parameter constraints at time instant T.
%
% See also: constraints/getStateConstraints,
% constraints/getAllConstraints, constraints/getInputConstraints, 
% constraints/getUnmeasurableInputConstraints,
% constraints/getReferenceConstraints

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

nx = object.nx;
nu = object.nu;
np = object.np;
nd = object.nd;
ny = object.ny;
nref = object.nref;

if nargin == 1
    
    H = cell(object.N+1,1);
    K = cell(object.N+1,1);
    
    for i = 1:object.N+1
        
        idx =  ismember(object.H{i}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
        H{i} = object.H{i}(idx,nx+nu+1:nx+nu+np);
        K{i} = object.K{i}(idx);
        
    end
    
    
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
        
        idx =  ismember(object.H{T+1}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
        
        H = object.H{T+1}(idx,nx+nu+1:nx+nu+np);
        K = object.K{T+1}(idx);
        
    elseif ischar(varargin{1})
        
        label = varargin{1};
        
        % Extracting only soft constraints
        if strcmpi(label,'soft')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx =  ismember(object.H{i}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
                idxsoft = object.soft{i} == 1;
                
                idx = idx & idxsoft;
                
                H{i} = object.H{i}(idx,nx+nu+1:nx+nu+np);
                K{i} = object.K{i}(idx);
                
            end
            
            % Extracting only hard constraints
        elseif strcmpi(label,'hard')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx =  ismember(object.H{i}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
                idxsoft = object.soft{i} == 0;
                
                idx = idx & idxsoft;
                
                H{i} = object.H{i}(idx,nx+nu+1:nx+nu+np);
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
        
        idx =  ismember(object.H{T+1}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
        idxsoft = object.soft{T+1} == 1;
        
        idx = idx & idxsoft;
        
        H = object.H{T+1}(idx,nx+nu+1:nx+nu+np);
        K = object.K{T+1}(idx);
        
        % Extracting only hard constraints at time T
    elseif strcmpi(label,'hard')
        
        idx =  ismember(object.H{T+1}(:,[1:nx+nu nx+nu+np+1:end]),zeros(1,nx+nu+nd+ny+nref),'rows');
        idxsoft = object.soft{T+1} == 0;
        
        idx = idx & idxsoft;
        
        H = object.H{T+1}(idx,nx+nu+1:nx+nu+np);
        K = object.K{T+1}(idx);
        
    else
        
        error('Label must be either ''soft'' or ''hard''');
        
    end
    
else
    
    error('Wrong input number');
end