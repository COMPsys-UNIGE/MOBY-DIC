function plot(object,f)
% disp   Plots the PWAG function
%
% plot(OBJ)
% Plots all the components of the vector PWAG function OBJ.
%
% plot(OBJ,IDX)
% Plots the components of the vector PWAG function specified in the
% vector of indices IDX.
%
% This function exploits MPT3 to perform the plots.

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

if ~exist('f','var')
    f = 1:object.ny;
end

if floor(f) ~= f
    error('f must be an integer number');
end
if any(f < 1) || any(f > object.ny)
    error(['f must be an integer number between 1 and ',num2str(object.ny)]);
end

if object.nx > 2
    error('Only functions defined over one- or two-dimensional domains can be plotted.');
end


figure
hold on
[Hd , Kd] = object.getDomain;

reg = object.getRegions;
P = [];
for i=1:numel(reg)
    Ptmp = Polyhedron(reg(i).H,reg(i).K);
    fun = AffFunction(reg(i).F,reg(i).G);
    Ptmp = Ptmp.addFunction(fun,'u');
    P = [P;Ptmp];
end

%P.plot();

% % Get regions
% regions = object.getRegions();
% 
% for i = 1:object.nr
%     H = regions(i).H;
%     K = regions(i).K;
%     F = regions(i).F(f,:);
%     G = regions(i).G(f);
%     
%     % Create Polyhedron object
%     P(i) = Polyhedron(H,K);
%     % Create affine function
%     fun = AffFunction(F, G);
%     % Add affine function to the polyhedron
%     P(i).addFunction(fun, 'fun');
% end

if object.nx == 1
    figure
    for i = 1:numel(f)
        subplot(1,numel(f),i)
        P.fplot('position',i);
        xlabel(object.xnames{1})
        ylabel(object.ynames{f(i)})
    end
end

if object.nx == 2
    figure
    for i = 1:numel(f)
        subplot(1,numel(f),i)
        P.fplot('position',i);
        xlabel(object.xnames{1})
        ylabel(object.xnames{2})
        zlabel(object.ynames{f(i)})
    end
end
