function optout = pwasCset(opts)
% pwasCset   Checks and sets the options for the generation of C code for
%            the pwasFunction object
%
% OPTOUT = pwasCset(OPTS)
% OPTS is the structure defining the options. See the documentation of 
% pwasFunction.generateC for a list of all structure fields.

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

if ~isfield(opts,'generateMain')
    opt.generateMain = true;
end

% Set or create folder in which to store the VHDL files
if isfield(opts,'folder')
    folder = opts.folder;
else
    folder = '';
end
if ~isempty(folder)
    if strcmp(folder(end),'\')
        folder(end) = '/';
    elseif ~strcmp(folder(end),'/')
        folder = [folder,'/'];
    end
else
    created = 0;
    i = 1;
    while ~created
        testfolder = [pwd,'/C_circuit_',num2str(i),'/'];
        if ~exist(testfolder,'dir')
            folder = testfolder;
            created = 1;
        else
            i = i+1;
        end
    end
end
optout.folder = folder;

if ~isfield(opts,'generateMain')
    optout.generateMain = 'true';
else
    optout.generateMain = opts.generateMain;
end