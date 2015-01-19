% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Plot circular scatterers on existing figure (for example on the scattered
% field). The figure is not cleared before drawing the disks.
% ----------------------------------------------
% PlotCirclesOnFigure(O, a, fig_index, ...)
% ----------------------------------------------
% Input arguments
% ---------------
%
% O         [2 x N_scat] : matrix of the centers of the disks, O(1,p)=x_coord of Op
% a         [1 x N_scat] : vector of the radii of the disks
% fig_index [1 x 1]      : index of the figure to draw circles on
%
% ----------------------------------------------
% Options
% -------
% PlotCirclesOnFigure(..., 'Color', COLOR)
% apply the COLOR color to lines (same as the plot function)
%
% PlotCirclesOnFigure(..., 'LineWidth', LINEWIDTH)
% set the line width to LINEWIDTH (same as the plot function)
%
% PlotCirclesOnFigure(..., 'zdata', ZDATA)
% set the zdata of the figure to ZDATA (same as the plot function)
%
% REMARK:
% -------
% When drawing disks on a figure containing some values, do not forget to
% set the zdata to the max value !
% examples: 
% > figure(1);
% > surf(X,Y,U);
% > PlotCirclesOnFigure(O, a, 1, 'zdata', max(max(U)))
%
% See also plot
%

function [] = PlotCircles(O, a, fig_index, varargin)
N_scat = length(a);
nvar = length(varargin);
thisColor = 'b';
thisLineWidth = ones(N_scat,1);
thisZdata = 100*ones(N_scat,1);
cpt = 1;
while (cpt < nvar)
    if(ischar(varargin{cpt}))
        if(strcmp(varargin{cpt}, 'LineWidth'))
            thisLineWidth = varargin{cpt+1};
            cpt = cpt + 2;
        elseif(strcmp(varargin{cpt}, 'zdata'))
            thisZdata =  varargin{cpt+1};
            cpt = cpt + 2;
        elseif(strcmp(varargin{cpt}, 'Color'))
            thisColor =  varargin{cpt+1};
            cpt = cpt + 2;
        else
            warning('Unknown option')
            cpt = cpt + 1;
        end
    end
end

if(isscalar(thisLineWidth))
   thisLineWidth = thisLineWidth*ones(N_scat,1); 
end
if(isscalar(thisZdata))
   thisZdata = thisZdata*ones(N_scat,1); 
end

figure(fig_index)
hold on
theta = [1:361]*2*pi/360;
nb = length(theta);
for p = 1:N_scat
    myzdata = thisZdata(p)*ones(1,nb);
    plot(O(1,p).*ones(1,nb) + a(p).*cos(theta), O(2,p).*ones(1,nb) + a(p).*sin(theta), thisColor, 'zdata', myzdata, 'LineWidth', thisLineWidth(p))
end
hold off;
