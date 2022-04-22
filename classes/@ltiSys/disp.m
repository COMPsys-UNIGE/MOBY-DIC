function disp(object)
% disp    Displays some information about the ltiSys object
%
% disp(OBJ)
% OBJ is the ltiSys object.

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


if object.nx == 0
    disp('Empty ltiSys object');
else
    if strcmp(object.mode,'ct')
        disp('Linear continuous-time system with: ')
        disp([' - ',num2str(object.nx),' state variables']);
        disp([' - ',num2str(object.nu),' input variables']);
        disp([' - ',num2str(object.ny),' output variables']);
        disp([' - ',num2str(object.np),' parameters']);
        disp([' - ',num2str(object.nd),' unmeasurable inputs']);
        disp(' ')
        
        lbx = repmat('[',object.nx,1);
        rbx = repmat(']',object.nx,1);
        lby = repmat('[',object.ny,1);
        rby = repmat(']',object.ny,1);
        xdot = repmat('        ',object.nx,1);
        xdot(floor(object.nx/2+1),:) = 'dx/dt = ';
        x = repmat('  ',object.nx,1);
        x(floor(object.nx/2+1),:) = ' x';
        xy = repmat('  ',object.ny,1);
        xy(floor(object.ny/2+1),:) = ' x';
        xys = repmat('  ',object.ny,1);
        xys(floor(object.ny/2+1),:) = ' x';
        ux = repmat('  ',object.nx,1);
        ux(floor(object.nx/2+1),:) = ' u';
        uy = repmat('  ',object.ny,1);
        uy(floor(object.ny/2+1),:) = ' u';
        px = repmat('  ',object.nx,1);
        px(floor(object.nx/2+1),:) = ' p';
        dx = repmat('  ',object.nx,1);
        dx(floor(object.nx/2+1),:) = ' d';
        py = repmat('  ',object.ny,1);
        py(floor(object.ny/2+1),:) = ' p';
        dy = repmat('  ',object.ny,1);
        dy(floor(object.ny/2+1),:) = ' d';
        y = repmat('        ',object.ny,1);
        y(floor(object.ny/2+1),:) = '    y = ';
        plusx = repmat('   ',object.nx,1);
        plusx(floor(object.nx/2+1),:) = ' + ';
        plusy = repmat('   ',object.ny,1);
        plusy(floor(object.ny/2+1),:) = ' + ';
        
        strx = [xdot lbx num2str(object.A) rbx x];
        stry = [y lby num2str(object.C) rby xy];
        
        if object.nu > 0
            strx = [strx plusx lbx num2str(object.B) rbx ux];
            stry = [stry plusy lby num2str(object.D) rby uy];
        end
        
        if object.np > 0
            strx = [strx plusx lbx num2str(object.Ex) rbx px];
            stry = [stry plusy lby num2str(object.Ey) rby py];
        end
        
        if object.nd > 0
            strx = [strx plusx lbx num2str(object.Fx) rbx dx];
            stry = [stry plusy lby num2str(object.Fy) rby dy];
        end
        
        strx = [strx plusx lbx num2str(object.Gx) rbx];
        stry = [stry plusy lby num2str(object.Gy) rby];
        
        disp(strx)
        disp(' ')
        disp(stry)
        disp(' ')
        
    else
        
        disp('Linear discrete-time system with: ')
        disp([' - ',num2str(object.nx),' state variables']);
        disp([' - ',num2str(object.nu),' input variables']);
        disp([' - ',num2str(object.ny),' output variables']);
        disp([' - ',num2str(object.np),' parameters']);
        disp([' - ',num2str(object.nd),' unmeasurable inputs']);
        disp(' ')
        disp(['Sampling time: ',num2str(object.Ts),' seconds']);
        disp(' ')
        
        lbx = repmat('[',object.nx,1);
        rbx = repmat(']',object.nx,1);
        lby = repmat('[',object.ny,1);
        rby = repmat(']',object.ny,1);
        xdot = repmat('         ',object.nx,1);
        xdot(floor(object.nx/2+1),:) = 'x(k+1) = ';
        x = repmat('     ',object.nx,1);
        x(floor(object.nx/2+1),:) = ' x(k)';
        xy = repmat('     ',object.ny,1);
        xy(floor(object.ny/2+1),:) = ' x(k)';
        xys = repmat('     ',object.ny,1);
        xys(floor(object.ny/2+1),:) = ' x(k)';
        ux = repmat('     ',object.nx,1);
        ux(floor(object.nx/2+1),:) = ' u(k)';
        uy = repmat('     ',object.ny,1);
        uy(floor(object.ny/2+1),:) = ' u(k)';
        px = repmat('  ',object.nx,1);
        px(floor(object.nx/2+1),:) = ' p';
        dx = repmat('  ',object.nx,1);
        dx(floor(object.nx/2+1),:) = ' d';
        py = repmat('  ',object.ny,1);
        py(floor(object.ny/2+1),:) = ' p';
        dy = repmat('  ',object.ny,1);
        dy(floor(object.ny/2+1),:) = ' d';
        y = repmat('        ',object.ny,1);
        y(floor(object.ny/2+1),:) = ' y(k) = ';
        plusx = repmat('   ',object.nx,1);
        plusx(floor(object.nx/2+1),:) = ' + ';
        plusy = repmat('   ',object.ny,1);
        plusy(floor(object.ny/2+1),:) = ' + ';
        
        strx = [xdot lbx num2str(object.A) rbx x];
        stry = [y lby num2str(object.C) rby xy];
        
        if object.nu > 0
            strx = [strx plusx lbx num2str(object.B) rbx ux];
            stry = [stry plusy lby num2str(object.D) rby uy];
        end
        
        if object.np > 0
            strx = [strx plusx lbx num2str(object.Ex) rbx px];
            stry = [stry plusy lby num2str(object.Ey) rby py];
        end
        
        if object.nd > 0
            strx = [strx plusx lbx num2str(object.Fx) rbx dx];
            stry = [stry plusy lby num2str(object.Fy) rby dy];
        end
        
        strx = [strx plusx lbx num2str(object.Gx) rbx];
        stry = [stry plusy lby num2str(object.Gy) rby];
        
        disp(strx)
        disp(' ')
        disp(stry)
        disp(' ')
        
    end
    
    xn = [];
    for i=1:object.nx-1
        xn = [xn object.xnames{i} ', '];
    end
    xn = [xn object.xnames{end}];
    disp(['State names: ',xn])
    
    if object.nu > 0
        un = [];
        for i=1:object.nu-1
            un = [un object.unames{i} ', '];
        end
        un = [un object.unames{end}];
        disp(['Input names: ',un])
    end
    
    yn = [];
    for i=1:object.ny-1
        yn = [yn object.ynames{i} ', '];
    end
    yn = [yn object.ynames{end}];
    disp(['Output names: ',yn])
    
    if object.np > 0
        pn = [];
        for i=1:object.np-1
            pn = [pn object.pnames{i} ', '];
        end
        pn = [pn object.pnames{end}];
        disp(['Parameter names: ',pn])
    end
    
    if object.nd > 0
        dn = [];
        for i=1:object.nd-1
            dn = [dn object.dnames{i} ', '];
        end
        dn = [dn object.dnames{end}];
        disp(['Unmeasurable input names: ',dn])
    end
end

if object.hasObserver
    disp(' ')
    disp('An observer is associated to this system');
end

if object.hasController
    disp(' ')
    disp('A controller is associated to this system');
end