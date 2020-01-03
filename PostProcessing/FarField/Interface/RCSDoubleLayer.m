% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the Radar Cross Section (RCS) of a double-layer potential, in the case of circular scatterers.
% The density needs to be decomposed on the Fourier bases
% -------------------------------------------------------------------------
% Double-Layer potential of density lambda :
% (mind the "minus" sign)
% u(x) = - \int_{Gamma} dn_y G(x,y) lambda(y) dy.
% Gamma = boundary of the scatterers.
% -------------------------------------------------------------------------
% Give the RCS of the double-layer potential:
% RCS = 10.*log10(2.*pi*(abs(Far_Field)).^2);
% -------------------------------------------------------------------------
%   RCS = RCSDoubleLayer(O, a, M_modes, k, Density, theta)
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat] : Vector of the centers of the scatterers
% a       [1 x N_scat] : Vector of the radius of the scatterers
% M_modes [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                        number of modes kept for the Fourier basis of the
%                        scatterer number "p".
% k           [1 x 1]  : Wavenumber in the vacuum
% Density     [1 x N]  : Coefficients in the Fourier bases of the
%                        density lambda.
%                        N = 2*sum(M_modes+1).
% theta [1 x N_angle]  : Angle of observation (rad).
%                        Note that theta could be a vector.
%
% Output Arguments :
% -----------------
% RCS   [1 x N_angle]  : Vector of the RCS
%
%
%%

function R = RCSDoubleLayer(O, a, M_modes, k, Density, theta)
    R = RCS(O, a, M_modes, k, Density, theta, [0,1]);
end