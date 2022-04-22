function plotTree(object)
% plotTree   Plots the binary search tree associated to the domain partition
%
% plotTree(OBJ)
% OBJ is the pwagFunction object.
%
% See also: pwagFunction/computeTree

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

if isempty(object.tree)
    disp('pwagFunction object does not have a tree')
else
    tree = object.tree;
    maxdepth = tree.maxDepth;
    meandepth = tree.meanDepth;
    balance = tree.balance;
    nnodes = tree.numberOfNodes;
    nodes = tree.nodes;
    
    ascissas=cell(maxdepth+1,1);
    
    % Inizialization
    for i=1:maxdepth+1
        ascissas{i}=zeros(1,2^(i-1));
    end
    
    % Generate points ascissas
    ascissas{maxdepth+1}=[0:2^maxdepth-1]/(2^maxdepth-1);
    
    for i=maxdepth:-1:1
        for j=1:2^(i-1)
            ascissas{i}(j)=0.5*(ascissas{i+1}(2*j-1)+ascissas{i+1}(2*j));
        end
    end
    
    info=struct('points',zeros(nnodes,2),'leaves',zeros(nnodes,1),'names',zeros(nnodes,1),'index',0,'prevName',0);
    
    info=treeStepPlot(nodes,nodes(1),info,ascissas,maxdepth);
    
    figure
    hold on
    title(['Number of nodes: ',num2str(nnodes),' - Mean depth: ',num2str(meandepth),' - Balance: ',num2str(balance)]);
    for i=1:nnodes
        if info.leaves(info.names(i)+1)==0
            plot(info.points(info.names(i)+1,1),info.points(info.names(i)+1,2),'ks','MarkerEdgeColor','k','MarkerFaceColor','r');
        else
            plot(info.points(info.names(i)+1,1),info.points(info.names(i)+1,2),'kd','MarkerEdgeColor','k','MarkerFaceColor','b');
        end
        
        if info.names(i)~=0
            if mod(info.names(i),2)==1
                plot([info.points(info.names(i)+1,1) info.points(0.5*(info.names(i)-1)+1,1)],...
                    [info.points(info.names(i)+1,2) info.points(0.5*(info.names(i)-1)+1,2)],'k');
            else
                plot([info.points(info.names(i)+1,1) info.points(0.5*(info.names(i)-2)+1,1)],...
                    [info.points(info.names(i)+1,2) info.points(0.5*(info.names(i)-2)+1,2)],'k');
            end
        end
    end
    axis tight
end
