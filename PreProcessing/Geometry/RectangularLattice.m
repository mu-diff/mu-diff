% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% TriangularLattice produces a set of center of obstacles (without radius) 
% periodically placed as a rectangular lattice. If "x" designes a center
% then the results looks like:
% x   x   x   x   x  
% x   x   x   x   x  
% x   x   x   x   x  
% x   x   x   x   x  
% x   x   x   x   x  
% x   x   x   x   x  
%
% The number of rows Ny and the number of obstacles Nx can be specified in 
% addition to the x-distance between two centers and the y-distance between
% two rows.
% By default, the first disks is located on (0,0) and is located on
% bottom-left of the lattice.
%
% -----------------------------------
% O = RectangularLattice(bx,by,Nx,Ny)
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
%   RectangularLattice(..., 'Origin', Ostart)
% place the first disk on Ostart position, Ostart being a [2x1] vector.
% Default: [0;0].
% 
%   RectangularLattice(..., 'Direction', +/- 1)
% place the row in the increasing y (+1) or decreasing y (-1). Default: +1
%
%   RectangularLattice(..., 'Center', Ocenter)
% center the collection on the point Ocenter.
% Default: none but erase Ostart ('Origin' option)
%
% See also TriangularLattice, CreateRandomDisks, PlotCircles
%

function O=RectangularLattice(bx, by, Nx, Ny, varargin)

% First point
Ostart=[0;0];
direction = 1; %from bottom to top

nvarargin = length(varargin);
cpt_arg = 1;
while(cpt_arg <= nvarargin)
   if(strcmp(varargin{cpt_arg}, 'Origin'))
       Ostart = varargin{cpt_arg+1};
       cpt_arg = cpt_arg +2;
   elseif(strcmp(varargin{cpt_arg}, 'Direction'))
       direction = varargin{cpt_arg+1};
       cpt_arg = cpt_arg +2;
   elseif(strcmp(varargin{cpt_arg}, 'Center'))
       Ocenter = varargin{cpt_arg+1};
       xOstart = Ocenter(1) - bx*(Nx-1)/2;
       yOstart = Ocenter(2) - by*(Ny-1)/2;
       Ostart = [xOstart; yOstart];
       cpt_arg = cpt_arg +2;
   else
       warning(['Unkown option ', varargin{cpt_arg}]);
       cpt_arg = cpt_arg +1;
   end
end

%Creation of the x-row of even or uneven index
x_vect=Ostart(1):bx:Ostart(1)+(Nx-1)*bx;

Ox =[];
Oy =[];

%Loop on every rows
for cpt = 0:Ny-1
    current_y = Ostart(2) + direction *by *cpt;
    Ox =[Ox, x_vect];
    Oy = [Oy, current_y*ones(1,Nx)];
end

O = [Ox;Oy];

