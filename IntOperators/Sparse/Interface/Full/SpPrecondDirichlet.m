% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the matrix of the integral operator obtained by preconditioning
% the single-layer operator by its single-scattering operator for a
% collection of sound-soft obstacles (Dirichlet), in the Fourier basis.
% SPARSE VERSION.
% -------------------------------------------------------------------------
%   A = SpPrecondDirichlet(O, a, M_modes, k)
%
% Output Arguments :
% -----------------
% A : cell(3,1) of the system
%
% Input Arguments (N_scat = length(a) = number of obstacles):
% -----------------
% O       [2 x N_scat] : Vector of the centers of the scatterers
% a       [1 x N_scat] : Vector of the radius of the scatterers
% M_modes [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                        number of modes kept for the Fourier basis of the
%                        scatterer number "p".
% k       [1 x 1]      : Wavenumber in the vacuum
%
% See also PrecondDirichlet, SpBlockIntegralOperator, SpIntegralOperator, SpMatVec

%%
function A = SpPrecondDirichlet(O, a, M_modes, k)
    A = SpIntegralOperator(O, a, M_modes, k, 6);
end
