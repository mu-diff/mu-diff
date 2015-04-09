% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the radar cross section of a solution computed with DirichletSolver, for the 
% Dirichlet problem of multiple scattering by disk.
% -------------------------------------------------------------------------
%   R = DirichletRCS(Solution, theta)
%
% OUTPUT ARGUMENTS:
% -----------------
% R  (sizeof(theta)) : Radar cross section of the scattered field
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by DirichletSolver
% theta      [1 x N]    : Angles of observation (in rad.)
%
% See also DirichletSolver, DirichletFarField, DirichletNearField
%
function R = DirichletRCS(Solution, theta)
F = DirichletFarField(Solution, theta);
R = FarField_to_RCS(F);
