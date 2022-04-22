function xupd = update(object,xcur,ucur,pcur,ycur)
% predict   Updates the system state based on current outputs
%
% XUPD = update(OBJ,XCUR,UCUR,PCUR,YCUR)
% XCUR is the estimated state at current time instant, UCUR, PCUR and YCUR 
% are the measured inputs, parameters and outputs at current time istant. 
% XUPD is the updated state for current time instant, based on the 
% measurements. OBJ is the kalmanFilter object.

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

% TO DO
% Controlla documentazione (anche di predict)


    if numel(xcur) ~= object.nx+object.nd
        error(['xcur must be a vector with ',num2str(object.nx+object.nd),' elements']);
    end
    if numel(ucur) ~= object.nu
        error(['u must be a vector with ',num2str(object.nu),' elements']);
    end
    if numel(pcur) ~= object.np
        error(['p must be a vector with ',num2str(object.np),' elements']);
    end
    if numel(ycur) ~= object.ny
        error(['y must be a vector with ',num2str(object.ny),' elements']);
    end

    % Get observer matrices
    [Aobs, Bobs, Cobs, Dobs, Gxobs, Gyobs] = object.getMatrices();
    
    % Take only last rows of matrices
    Cobs = Cobs(object.ny+1:end,:);
    Dobs = Dobs(object.ny+1:end,:);
    Gyobs = Gyobs(object.ny+1:end,:);
    
    % Compute predicted state
    xupd = Cobs*xcur(:)+Dobs*[ucur(:);pcur(:);ycur(:)]+Gyobs;
    

end