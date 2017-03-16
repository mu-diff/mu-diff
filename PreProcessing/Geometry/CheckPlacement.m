% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Check if the list of obstacles sent are well placed with respect to parameters.
%res = CheckPlacement(O, a) : will check if the list of disks of centers O
%                             and radii a do not overlap (=true)
% Input arguments:
% ----------------
%   
% O  [2, N_scat]    : Coordinates of the N_scat disks
% a  [1, N_scat]    : Radii of the circular obstacles
%
% Output arguments:
% -----------------
% res  [1x1]        : true if obstacles are well placed, false otherwise
%
% Options 1: res = CheckPlacement(O, a, d_min, d_max)
% ---------------------------------------------------
%   Will check if two disks are distant from at least d_min and at most
%   from d_max
%   dmin(=realmin)  [1, 1]    : minimal authorized distance between two disks
%   dmax(=realmax)  [1, 1]    : maximal authorized distance between two disks
%
% Options 2: res = CheckPlacement(O, a, d_min, d_max, Onew, anew)
% ---------------------------------------------------------------
%   Will check if the disk (Onew, anew) is well placed compare to all other
%   disks, with the good distances. Note that other disks are then NOT
%   checked.
%   Onew  [2, 1]    : Coordinate of the new disk
%   anew  [1, 1]    : Radius of the new disk
%
% Remarks:
% --------
% - If Onew is specified without anew, then Onew is not considered
% - If Onew and anew are specified then the disks (O,a) are NOT checked.
%   They are assumed to be well placed.
% - Be carreful with dmax. It is not the distance between the closest disk,
%   but with ALL OTHER disks.


function res = CheckPlacement(O, a, varargin)

res = true;
N_scat = length(a);
if(size(O,2) ~= N_scat)
    error('O and a should have the same number of columns');
end

CheckAllObstacles = true;

dmin = realmin;
dmax = realmax;

if(nargin >= 3)
    dmin = varargin{1};
end
if(nargin >= 4)
    dmax = varargin{2};
end
if(nargin >= 6)
    Onew = varargin{3};
    anew = varargin{4};
    CheckAllObstacles = false;
end

if(dmin <= 0)
    dmin = realmin;
end
if(dmax <= 0)
    dmax = realmax;
end


if(CheckAllObstacles) %check all obstacles of center O and radii a
    for p =2:N_scat
        if(res>0)
           break;
        end
        Op = O(:,p);
        ap = a(p);
        res_min = any(sqrt((O(1,1:p-1) - Op(1)).^2+ (O(2,1:p-1) - Op(2)).^2) - a(1:p-1) - ap < dmin);
        res_max = any(sqrt((O(1,1:p-1) - Op(1)).^2+ (O(2,1:p-1) - Op(2)).^2) - a(1:p-1) - ap > dmax);
        res = res_min+res_max;
    end
else % check only new obstacle with obstacles of center 0 and radii a
    res_min = any(sqrt((O(1,:) - Onew(1)).^2+ (O(2,:) - Onew(2)).^2) - a - anew < dmin);
    res_max = any(sqrt((O(1,:) - Onew(1)).^2+ (O(2,:) - Onew(2)).^2) - a - anew > dmax);
    res = res_min+res_max;
end
res = ~res;