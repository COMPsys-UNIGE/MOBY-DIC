function [H, K] = getReferenceConstraints(object,varargin)
% getReferenceConstraints   Gets the H and K matrices defining reference 
%                           state or output constraints
%
% [H, K] = getReferenceConstraints(OBJ)
% Gets all constraints involving system reference states or outputs. 
% H and K are cell arrays, whose elements corresponds to different time 
% instants (up to time horizon N). Each elements of the cell arays defines 
% a system of inequalities in the form:
%
% H ref  <= K
%
% being ref the system reference states or outputs. OBJ is the constraints 
% object.
%
% [H, K] = getReferenceConstraints(OBJ,T)
% Gets all constraints involving system reference states or outputs at time 
% instant T. T must be a scalar between 0 and the time horizon N. 
% H and K are matrices defining all reference state or output constraints 
% in the form:
%
% H ref  <= K
%
% [H, K] = getReferenceConstraints(OBJ,'soft')
% [H, K] = getReferenceConstraints(OBJ,'hard')
% Gets only soft or hard reference state or output constraints.
%
% [H, K] = getReferenceConstraints(OBJ,'soft',T)
% [H, K] = getReferenceConstraints(OBJ,'hard',T)
% Gets only soft or hard reference state or output constraints at time 
% instant T.
%
% See also: constraints/getStateConstraints,
% constraints/getAllConstraints, constraints/getParameterConstraints, 
% constraints/getInputConstraints, constraints/getUnmeasurableInputConstraints

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

if nargin == 1
    
    H = cell(object.N+1,1);
    K = cell(object.N+1,1);
    
    for i = 1:object.N+1
        
        idx =  ismember(object.H{i}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
        
        H{i} = object.H{i}(idx,nx+nu+np+nd+ny+1:end);
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
        
        idx =  ismember(object.H{T+1}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
        
        H = object.H{T+1}(idx,nx+nu+np+nd+ny+1:end);
        K = object.K{T+1}(idx);
        
    elseif ischar(varargin{1})
        
        label = varargin{1};
        
        % Extracting only soft constraints
        if strcmpi(label,'soft')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx =  ismember(object.H{i}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
                idxsoft = object.soft{i} == 1;
                
                idx = idx & idxsoft;
                
                H{i} = object.H{i}(idx,nx+nu+np+nd+ny+1:end);
                K{i} = object.K{i}(idx);
                
            end
            
            % Extracting only hard constraints
        elseif strcmpi(label,'hard')
            
            H = cell(object.N+1,1);
            K = cell(object.N+1,1);
            for i = 1:object.N+1
                
                idx =  ismember(object.H{i}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
                idxsoft = object.soft{i} == 0;
                
                idx = idx & idxsoft;
                
                H{i} = object.H{i}(idx,nx+nu+np+nd+ny+1:end);
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
        
        idx =  ismember(object.H{T+1}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
        idxsoft = object.soft{T+1} == 1;
        
        idx = idx & idxsoft;
        
        H = object.H{T+1}(idx,nx+nu+np+nd+ny+1:end);
        K = object.K{T+1}(idx);
        
        % Extracting only hard constraints at time T
    elseif strcmpi(label,'hard')
        
        idx =  ismember(object.H{T+1}(:,1:nx+nu+np+nd+ny),zeros(1,nx+nu+np+nd+ny),'rows');
        idxsoft = object.soft{T+1} == 0;
        
        idx = idx & idxsoft;
        
        H = object.H{T+1}(idx,nx+nu+np+nd+ny+1:end);
        K = object.K{T+1}(idx);
        
    else
        
        error('Label must be either ''soft'' or ''hard''');
        
    end
    
else
    
    error('Wrong input number');
end