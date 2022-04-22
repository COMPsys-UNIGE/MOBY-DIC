function object = setMatrices(object,dyn,varargin)
% setMatrices    Sets the matrices defining the LTI system
%
% OBJ = setMatrices(OBJ,DYN,NAME,MATR)
%
% Sets the matrices defining the dynamics indexed by DYN of the PWA system. 
% NAME is a string identifying the matrix (possible values: 'A', 'B', 'C', 
% 'D', 'Ex', 'Ey', 'Fx', 'Fy', 'Gx' or 'Gy') and MATR is the corresponding 
% matrix. OBJ is the pwaSys object.

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



if nargin ~= 4
    error('Wrong number of inputs')
end

if dyn > object.nDyn
    error(['Dynamics ',num2str(dyn),' not created yet!']);
end

name = varargin{1};
matr = varargin{2};

if strcmpi(name,'A')
    A = matr;
    if size(A,1) ~= size(A,2)
        error('A matrix must be square');
    end
    if size(A,1) ~= object.nx
        error('A matrix is not compatible with number of states');
    end
    
    object.A{dyn} = A;
       
elseif strcmpi(name,'B')
    B = matr;
    if size(B,1) ~= object.nx
        error('B matrix is not compatible with number of states');
    end
    if size(B,2) ~= object.nu
        error('B matrix is not compatible with number of inputs');
    end
    
    object.B{dyn} = B;
    
elseif strcmpi(name,'C')
    C = matr;
    if size(C,1) ~= object.ny
        error('C matrix is not compatible with number of outputs');
    end
    if size(C,2) ~= object.nx
        error('C matrix is not compatible with number of states');
    end
    
    object.C{dyn} = C;
    
elseif strcmpi(name,'D')
    D = matr;
    if size(D,1) ~= object.ny
        error('D matrix is not compatible with number of outputs');
    end
    if size(D,2) ~= object.nu
        error('D matrix is not compatible with number of inputs');
    end
    
    object.D{dyn} = D;    
    
elseif strcmpi(name,'Ex')
    Ex = matr;
    if size(Ex,1) ~= object.nx
        error('E_x matrix is not compatible with number of states');
    end
    if size(Ex,2) ~= object.np
        error('E_x matrix is not compatible with number of parameters');
    end
    
    object.Ex{dyn} = Ex;   
    
elseif strcmpi(name,'Ey')
    Ey = matr;
    if size(Ey,1) ~= object.ny
        error('E_y matrix is not compatible with number of outputs');
    end
    if size(Ey,2) ~= object.np
        error('E_y matrix is not compatible with number of parameters');
    end
    
    object.Ey{dyn} = Ey;       
    
elseif strcmpi(name,'Fx')
    Fx = matr;
    if size(Fx,1) ~= object.nx
        error('F_x matrix is not compatible with number of states');
    end
    if size(Fx,2) ~= object.nd
        error('F_x matrix is not compatible with number of unmeasurable inputs');
    end
    
    object.Fx{dyn} = Fx;   
    
elseif strcmpi(name,'Fy')
    Fy = matr;
    if size(Fy,1) ~= object.ny
        error('F_y matrix is not compatible with number of outputs');
    end
    if size(Fy,2) ~= object.nd
        error('F_y matrix is not compatible with number of unmeasurable inputs');
    end
    
    object.Fy{dyn} = Fy;  
    
elseif strcmpi(name,'Gx')
    Gx = matr;
    if size(Gx,1) ~= object.nx
        error('G_x matrix is not compatible with number of states');
    end
    if size(Gx,2) ~= 1
        error('G_x must be a column vector');
    end
    
    object.Gx{dyn} = Gx;   
    
elseif strcmpi(name,'Gy')
    Gy = matr;
    if size(Gy,1) ~= object.ny
        error('G_y matrix is not compatible with number of outputs');
    end
    if size(Gy,2) ~= 1
        error('G_x must be a column vector');
    end
    
    object.Gy{dyn} = Gy;  
    
else 
    error('Wrong value for ''name''');
    
end
