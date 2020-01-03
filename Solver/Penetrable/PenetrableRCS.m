% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the radar cross section of a solution computed with PenetrableSolver, for the 
% penetrable problem of multiple scattering by disk.
% -------------------------------------------------------------------------
%   R = PenetrableRCS(Solution, theta)
%
% OUTPUT ARGUMENTS:
% -----------------
% R  (sizeof(theta)) : Radar cross section of the scattered field
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by PenetrableSolver
% theta      [1 x N]    : Angles of observation (in rad.)
%
% See also PenetrableSolver, PenetrableFarField, PenetrableNearField
%
function R = PenetrableRCS(Solution, theta)
F = PenetrableFarField(Solution, theta);
R = FarField_to_RCS(F);
