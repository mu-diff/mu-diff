% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Remove some centers of disks from the matrix O
% ----------------------------------------------------------
% [O,a] = RemoveDisk(O, ...)
% 
% Input (N_scat = number of obstacle = length(a))
% -----
%   O_old [2 x N_scat_old]  : Matrix of the coordinates of the center of the
%                             scatterers
%   a_old [1 x N_scat_old]  : Vector of the radii of the disks
%
% Output
% ------
%   O [2 x N_scat]  : Matrix of the coordinates of the center of the
%                     scatterers
%   a [1 x N_scat_old]  : Vector of the radii of the disks
%
% Options:
% --------
%
% [O,a] = RemoveDisk(..., 'X', [X1, X2, ..., XN])
% Remove all the points with X abscissa X1, X2, ..., or XN
%
% [O,a] = RemoveDisk(..., 'Y', [Y1, Y2, ..., YN])
% Remove all the points with Y ordinate Y1, Y2, ..., or YN
%
% [O,a] = RemoveDisk(..., 'XY', [[X1;Y1], [X2;Y2], ..., [XN;YN]])
% Remove all the points [X1;Y1], [X2;Y2], ..., and [XN;YN]
%
% [O,a] = RemoveDisk(..., 'Radius', [a1, a2, ..., aN])
% Remove all the disk with radius a1, a2, ..., or aN
%
% [O,a] = RemoveDisk(..., 'Verbosity', VERBOSITY)
% set VERBOSITY to 0 to avoid display message, to 1 to only show results,
% and to > 1 to see everything (default).
%
% See also CreateRandomDisks, RectangularLattice, TriangularLattice

function [O,a] = RemoveDisk(O_old, a_old, varargin)

Tol = 10^(-10);

X_to_avoid = [];
Y_to_avoid = [];
XY_to_avoid = [];
a_to_avoid = [];
nvarargin = length(varargin);
cpt_arg = 1;
VERBOSITY =2;
while(cpt_arg <= nvarargin)
   if(strcmp(varargin{cpt_arg}, 'X'))
    X_to_avoid = [X_to_avoid, varargin{cpt_arg+1}];
    cpt_arg = cpt_arg+2;
   elseif(strcmp(varargin{cpt_arg}, 'Y'))
    Y_to_avoid = [Y_to_avoid, varargin{cpt_arg+1}];
    cpt_arg = cpt_arg+2;
   elseif(strcmp(varargin{cpt_arg}, 'XY'))
    XY = varargin{cpt_arg+1};
    XY_to_avoid = [XY_to_avoid, [XY(1);XY(2)]];
    cpt_arg = cpt_arg+2;
   elseif(strcmp(varargin{cpt_arg}, 'Radius'))
    a_to_avoid = [a_to_avoid, varargin{cpt_arg+1}];
    cpt_arg = cpt_arg+2;
   elseif(strcmp(varargin{cpt_arg}, 'Verbosity'))
    VERBOSITY = varargin{cpt_arg+1};
    cpt_arg = cpt_arg+2;
   end
end

N_scat_old = size(O_old,2);
N_scat_delete = 0;
N_scat = 0;
O = [];
a = [];

for p=1:N_scat_old
    Xold = O_old(1,p);
    Yold = O_old(2,p);
    aold = a_old(p);
    isOK = true;
    for cpt=1:size(X_to_avoid,2)
       if(abs(Xold - X_to_avoid(cpt))<Tol)
           isOK = false;
           break;
       end
    end
    if(isOK)
        for cpt=1:size(Y_to_avoid,2)
           if(abs(Yold-Y_to_avoid(cpt))<Tol)
               isOK = false;
               break;
           end
        end
    end
    if(isOK)
        for cpt=1:size(XY_to_avoid,2)
           if(abs(Xold-XY_to_avoid(1, cpt))<Tol && abs(Yold-XY_to_avoid(2, cpt))<Tol)
               isOK = false;
               break;
           end
        end
    end
    if(isOK)
        for cpt=1:size(a_to_avoid,2)
           if(abs(aold - a_to_avoid(p))<Tol)
               isOK = false;
               break;
           end
        end
    end
    if(isOK)
       O = [O, [Xold;Yold]]; 
       a = [a, aold];
       N_scat = N_scat +1 ;
    else
       if(VERBOSITY > 1)
            disp(['Removing disk [',num2str(Xold),',',num2str(Yold),']']);
       end
       N_scat_delete = N_scat_delete+1;
    end
end

if(VERBOSITY > 0)
    disp(['Removed disks: ',num2str(N_scat_delete)]);
    disp(['Old number of disks: ',num2str(N_scat_old)]);
    disp(['New number of disks: ',num2str(N_scat)]);
end

end