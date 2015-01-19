% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% TriangularLattice produces a set of center of obstacles (without radius) 
% periodically placed as a triangular lattice. If "x" designes a center
% then the results looks like:
% x   x   x   x   x  
%   x   x   x   x   x
% x   x   x   x   x  
%   x   x   x   x   x
% x   x   x   x   x  
%   x   x   x   x   x
%
% The number of rows Ny and the number of obstacles Nx can be specified in 
% addition to the x-distance between two centers and the y-distance between
% two rows.
% By default, the first disks is located on (0,0) and is located on
% bottom-left of the lattice.
%
% -----------------------------------
% O = TriangularLattice(bx,by,Nx,Ny)
%
% OUTPUT ARGUMENTS:
% -----------------
% O [2 x N_scat]    : Matrix of the centers of the disks
%
% INPUT ARGUMENTS:
% ----------------
% bx [1x1]  : x-distance between two centers
% by [1x1]  : y-distance between two row of obstacles
% Nx [1x1]  : Number of obstacles on a row
% Ny [1x1]  : Number of rows
%
% OPTIONS:
% --------
%   TriangularLattice(..., 'Center', Ostart)
% place the first disk on Ostart position, Ostart being a [2x1] vector.
% Default: [0;0].
% 
%   TriangularLattice(..., 'Direction', +/- 1)
% place the row in the increasing y (+1) or decreasing y (-1). Default: +1
%
% See also RectangularLattice, CreateRandomDisks, PlotCircles
%

function O=TriangularLattice(bx, by, Nx, Ny, varargin)

% First point
Ostart=[0;0];
direction = 1; %from bottom to top

nvarargin = length(varargin);
cpt_arg = 1;
while(cpt_arg <= nvarargin)
   if(strcmp(varargin{cpt_arg}, 'Center'))
       Ostart = varargin{cpt_arg+1};
       cpt_arg = cpt_arg +2;
   elseif(strcmp(varargin{cpt_arg}, 'Direction'))
       direction = varargin{cpt_arg+1};
       cpt_arg = cpt_arg +2;
   else
       warning('Unkown option');
       cpt_arg = cpt_arg +1;
   end
end

%Creation of the x-row of even or uneven index
x_even=Ostart(1):bx:Ostart(1)+(Nx-1)*bx;
x_uneven=Ostart(1)+bx/2:bx:Ostart(1)+bx/2+(Nx-2)*bx;

Ox =[];
Oy =[];

%Loop on every rows
for cpt = 0:Ny-1
    current_y = Ostart(2) + direction *by *cpt;
    if(mod(cpt,2)==0)
        Ox =[Ox, x_even];
        Oy = [Oy, current_y*ones(1,Nx)];
    else
        Ox =[Ox, x_uneven];
        Oy = [Oy, current_y*ones(1,Nx-1)];
    end
end

O = [Ox;Oy];
