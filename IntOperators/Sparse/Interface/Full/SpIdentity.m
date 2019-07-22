% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the matrix of the identity operator for the multiple scattering
% problem in the Fourier bases.
% SPARSE VERSION.
% -------------------------------------------------------------------------
%   I = SpIdentity(O, a, M_modes)
%
% Output Arguments :
% -----------------
% I : cell(3,1) of the system
%
% Input Arguments (N_scat = length(a) = number of obstacles):
% -----------------
% O       [2 x N_scat] : Vector of the centers of the scatterers
% a       [1 x N_scat] : Vector of the radius of the scatterers
% M_modes [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                        number of modes kept for the Fourier basis of the
%                        scatterer number "p".
%
% See also SpBlockIntegralOperator, SpIntegralOperator, SpMatVec

%%
function I = SpIdentity(O, a, M_modes)
    I = SpIntegralOperator(O, a, M_modes, 1, 1);
end