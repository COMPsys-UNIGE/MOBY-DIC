function object = setConstraints(object,varargin)
% setConstraints   Sets inequality hard constraints
%
% OBJ = setConstraints(OBJ,LABEL,LBOUND,UBOUND)
% Sets hard saturation constraints in the form:
%
%   LBOUND <= X <= UBOUND
%
% LABEL can assume values 'x', 'u', 'y', 'p', 'd', 'r' followed by a number
% (e.g. 'x1', 'u2', 'y1', 'p3', 'd2', 'r1').
% - If LABEL is 'xj', state saturation constraints in the form
%   LBOUND <= xj <= UBOUND are set for j-th state variables.
% - If LABEL is 'uj', saturation constraints in the form
%   LBOUND <= uj <= UBOUND are set for j-th input variables.
% - If LABEL is 'yj', saturation constraints in the form
%   LBOUND <= yj <= UBOUND are set for j-th output variables.
% - If LABEL is 'pj', saturation constraints in the form
%   LBOUND <= pj <= UBOUND are set for j-th parameter.
% - If LABEL is 'dj', saturation constraints in the form
%   LBOUND <= dj <= UBOUND are set for j-th unmeasurable input.
% - If LABEL is 'xrefj', saturation constraints in the form
%   LBOUND <= rj <= UBOUND are set for j-th reference state.
% - If LABEL is not followed by a number and LBOUND and UBOUND are
%   vectors, the constraints are imposed to all components of states, inputs,
%   outputs, parameters, unmeasurable inputs or reference states.
% The constraints are imposed at all time instants. In particular,
% constraints involving only parameters, unmeasurable inputs or reference
% states are imposed only at time k, constraints involving also states ans/or
% outputs are imposed from time k to time k+N, constraints involving also
% inputs are imposed from time k to time k+N-1.
%
% OBJ = setConstraints(OBJ,LABEL,LBOUND,UBOUND,T)
% As before, but the constraints are imposed at time instant T.
%
% OBJ = setConstraints(OBJ,H,K)
% Sets hard mixed constraints in the form
%    _   _
%   |  x   |
%   |  u   |
% H |  p   | <= K
%   |  d   |
%   |  y   |
%   |_xref_|
%
% The constraints are imposed at all time instants, up to time horizon N.
%
% OBJ = setConstraints(OBJ,H,K,T)
% As before, but the constraints are imposed at time instant T.
%
% See also: constraints/setSoftConstraints.

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
N = object.N;

if nargin == 3
    label = 'm';
    H = varargin{1};
    K = varargin{2};
    T = [];
elseif nargin == 4
    if ischar(varargin{1})
        label = varargin{1};
        lbound = varargin{2};
        ubound = varargin{3};
        % Transform in column vectors
        lbound = lbound(:);
        ubound = ubound(:);
        % Check if ubound us greater than lbound
        if any(lbound > ubound)
            error('Lower bound is greater than upper bound!');
        end
        T = [];
    else
        label = 'm';
        H = varargin{1};
        K = varargin{2};
        T = varargin{3};
    end
elseif nargin == 5
    label = varargin{1};
    lbound = varargin{2};
    ubound = varargin{3};
    % Transform in column vectors
    lbound = lbound(:);
    ubound = ubound(:);
    % Check if ubound us greater than lbound
    if any(lbound > ubound)
        error('Lower bound is greater than upper bound!');
    end
    T = varargin{4};
else
    error('Wrong input arguments');
end

% This happens if label = 'x', 'u', 'y', 'p', 'd', 'r' or 'm'
if numel(label) == 1
    
    % State constraints
    if strcmpi(label(1),'x')
        
        % State saturation constraints
        if numel(lbound) == nx && numel(ubound) == nx
            Hx = [-eye(nx);eye(nx)];
            Kx = [-lbound;ubound];
            Hu = zeros(size(Hx,1),nu);
            Hy = zeros(size(Hx,1),ny);
            Hp = zeros(size(Hx,1),np);
            Hd = zeros(size(Hx,1),nd);
            Hr = zeros(size(Hx,1),nref);
        else
            error(['lbound and ubound must be vectors with ', num2str(nx),' elements']);
        end
        
        if isempty(T)
            T = 0:N;
        end
        
        if any(T < 0) || any(T > N)
            error('Wrong value for T. State constraints can be imposed from instant k to k+N');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Kx];
            
            % At instant k, hard constraints are imposed (to define the
            % domain)
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
            
        end
        
        % Input constraints
    elseif strcmpi(label(1),'u')
        
        % State saturation constraints
        if numel(lbound) == nu && numel(ubound) == nu
            Hu = [-eye(nu);eye(nu)];
            Ku = [-lbound;ubound];
            Hx = zeros(size(Hu,1),nx);
            Hy = zeros(size(Hu,1),ny);
            Hp = zeros(size(Hu,1),np);
            Hd = zeros(size(Hu,1),nd);
            Hr = zeros(size(Hu,1),nref);
        else
            error(['lbound and ubound must be vectors with ', num2str(nu),' elements']);
        end
        
        if isempty(T)
            T = 0:N-1;
        end
        
        if any(T < 0) || any(T > N-1)
            error('Wrong value for T. Input constraints can be imposed from instant k to k+N-1');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Ku];
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
        end
        
        % Parameter constraints
    elseif strcmpi(label(1),'p')
        
        % State saturation constraints
        if numel(lbound) == np && numel(ubound) == np
            Hp = [-eye(np);eye(np)];
            Kp = [-lbound;ubound];
            Hx = zeros(size(Hp,1),nx);
            Hu = zeros(size(Hp,1),nu);
            Hy = zeros(size(Hp,1),ny);
            Hd = zeros(size(Hp,1),nd);
            Hr = zeros(size(Hp,1),nref);
        else
            error(['lbound and ubound must be vectors with ', num2str(np),' elements']);
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Parameter constraints can be imposed only on instant k');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Kp];
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
        end
        
        % Unmeasurable input constraints
    elseif strcmpi(label(1),'d')
        
        % State saturation constraints
        if numel(lbound) == nd && numel(ubound) == nd
            Hd = [-eye(nd);eye(nd)];
            Kd = [-lbound;ubound];
            Hx = zeros(size(Hd,1),nx);
            Hy = zeros(size(Hd,1),ny);
            Hu = zeros(size(Hd,1),nu);
            Hp = zeros(size(Hd,1),np);
            Hr = zeros(size(Hd,1),nref);
        else
            error(['lbound and ubound must be vectors with ', num2str(nd),' elements']);
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Unmeasurable input constraints can be imposed only on instant k');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Kd];
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
        end
        
        % Output constraints
    elseif strcmpi(label(1),'y')
        
        % Output saturation constraints
        if numel(lbound) == ny && numel(ubound) == ny
            Hy = [-eye(ny);eye(ny)];
            Ky = [-lbound;ubound];
            Hx = zeros(size(Hy,1),nx);
            Hd = zeros(size(Hy,1),nd);
            Hu = zeros(size(Hy,1),nu);
            Hp = zeros(size(Hy,1),np);
            Hr = zeros(size(Hy,1),nref);
        else
            error(['lbound and ubound must be vectors with ', num2str(ny),' elements']);
        end
        
        if isempty(T)
            T = 0:N;
        end
        
        if any(T < 0) || any(T > N)
            error('Wrong value for T. Output constraints can be imposed from instant k to k+N');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Ky];
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
        end
        
    elseif strcmpi(label(1),'r')
        
        % State saturation constraints
        if numel(lbound) == nref && numel(ubound) == nref
            Hr = [-eye(nref);eye(nref)];
            Kr = [-lbound;ubound];
            Hx = zeros(size(Hr,1),nx);
            Hu = zeros(size(Hr,1),nu);
            Hy = zeros(size(Hr,1),ny);
            Hp = zeros(size(Hr,1),np);
            Hd = zeros(size(Hr,1),nd);
        else
            error(['lbound and ubound must be vectors with ', num2str(nref),' elements']);
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Reference state constraints can be imposed only on instant k');
        end
        
        for i = T+1
            object.H{i} = [object.H{i}; Hx Hu Hp Hd Hy Hr];
            object.K{i} = [object.K{i}; Kr];
            object.soft{i} = [object.soft{i};zeros(size(Hx,1),1)];
        end
        
        % Mixed constraints
    elseif strcmpi(label(1),'m')
        
        if size(H,2) == nx+nu+np+nd+ny+nref && size(H,1) == size(K,1) && size(K,2) == 1
            
            Hu = H(:,nx+1:nx+nu);
            
            % If mixed constraints also involve u
            if any(Hu)
                
                if isempty(T)
                    T = 0:N-1;
                end
                
                if any(T < 0) || any(T > N-1)
                    error('Wrong value for T. Constraints involving inputs can be imposed from instant k to k+N-1');
                end
                
                for i = T+1
                    object.H{i} = [object.H{i}; H];
                    object.K{i} = [object.K{i}; K];
                    object.soft{i} = [object.soft{i};zeros(size(H,1),1)];
                end
                
            else
                
                if isempty(T)
                    T = 0:N;
                end
                
                if any(T < 0) || any(T > N)
                    error('Wrong value for T. Constraints can be imposed from instant k to k+N');
                end
                
                for i = T+1
                    object.H{i} = [object.H{i}; H];
                    object.K{i} = [object.K{i}; K];
                    object.soft{i} = [object.soft{i};zeros(size(H,1),1)];
                end
                
            end
            
        else
            error('Wrong dimensions for matrices H and K.');
        end
        
    else
        error('Wrong value for ''label''');
    end
    
    
    % This happens if label is followed by a number (e.g. x1, u2, p4)
else
    
    num = str2double(label(2:end));
    
    % State constraints
    if strcmpi(label(1),'x')
        
        if num < 1 || num > nx
            error('State not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hx = [-1;1];
            Kx = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0:N;
        end
        
        if any(T < 0) || any(T > N)
            error('Wrong value for T. State constraints can be imposed from instant k to k+N');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,num) = Hx;
        Ktmp = Kx;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
        % Input constraints
    elseif strcmpi(label(1),'u')
        
        if num < 1 || num > nu
            error('Input not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hu = [-1;1];
            Ku = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0:N-1;
        end
        
        if any(T < 0) || any(T > N-1)
            error('Wrong value for T. Input constraints can be imposed from instant k to k+N-1');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,nx+num) = Hu;
        Ktmp = Ku;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
        % Parameter constraints
    elseif strcmpi(label(1),'p')
        
        if num < 1 || num > np
            error('Parameter not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hp = [-1;1];
            Kp = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Parameters constraints can only be imposed at instant k');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,nx+nu+num) = Hp;
        Ktmp = Kp;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
        % Unmeasurable input constraints
    elseif strcmpi(label(1),'d')
        
        if num < 1 || num > nd
            error('Unmeasurable input not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hd = [-1;1];
            Kd = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Unmeasurable input constraints can only be imposed at instant k');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,nx+nu+np+num) = Hd;
        Ktmp = Kd;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
        % Output constraints
    elseif strcmpi(label(1),'y')
        
        if num < 1 || num > ny
            error('Output not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hy = [-1;1];
            Ky = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0:N;
        end
        
        if any(T < 0) || any(T > N)
            error('Wrong value for T. Output constraints can be imposed from instant k to k+N');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,nx+nu+np+nd+num) = Hy;
        Ktmp = Ky;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
        
    elseif strcmpi(label(1),'r')
        
        if num < 1 || num > nref
            error('Reference state not available');
        end
        
        if numel(lbound) == 1 && numel(ubound) == 1
            
            Hr = [-1;1];
            Kr = [-lbound;ubound];
            
        else
            error('lbound and ubound must be scalars');
        end
        
        if isempty(T)
            T = 0;
        end
        
        if any(T ~= 0)
            error('Wrong value for T. Unmeasurable input constraints can only be imposed at instant k');
        end
        
        Htmp = zeros(2,nx+nu+np+nd+ny+nref);
        Htmp(:,nx+nu+np+nd+ny+num) = Hr;
        Ktmp = Kr;
        
        for i = T+1
            object.H{i} = [object.H{i}; Htmp];
            object.K{i} = [object.K{i}; Ktmp];
            object.soft{i} = [object.soft{i};zeros(size(Htmp,1),1)];
        end
        
    else
        
        error('Wrong value for ''label''');
    end
    
end


