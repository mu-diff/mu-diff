% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the radar cross section of a solution computed with NeumannSolver, for the 
% Neumann problem of multiple scattering by disk.
% -------------------------------------------------------------------------
%   R = NeumannRCS(Solution, theta)
%
% OUTPUT ARGUMENTS:
% -----------------
% R  (sizeof(theta)) : Radar cross section of the scattered field
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by NeumannSolver
% theta      [1 x N]    : Angles of observation (in rad.)
%
% See also NeumannSolver, NeumannFarField, NeumannNearField
%
function R = NeumannRCS(Solution, theta)
F = NeumannFarField(Solution, theta);
R = FarField_to_RCS(F);
