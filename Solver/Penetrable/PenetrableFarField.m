% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the far field of a solution computed with PenetrableSolver, for the 
% penetrable problem of multiple scattering by disk.
% -------------------------------------------------------------------------
%   F = PenetrableFarField(Solution, theta)
%
% OUTPUT ARGUMENTS:
% -----------------
% F  (sizeof(theta)) : Far field of the scattered field
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by PenetrableSolver
% theta      [1 x N]    : Angles of observation (in rad.)
%
% See also PenetrableSolver, PenetrableRCS, PenetrableNearField
%
function F = PenetrableFarField(Solution, theta)
O = Solution{1};
a = Solution{2};
k = Solution{3};
M_modes = Solution{5};
rho_plus = Solution{8};
%% Far Field
F = FarField(O, a, M_modes, k, theta, rho_plus, [1,0]);
end