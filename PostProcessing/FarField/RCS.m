% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the Radar Cross Section of a combination of single- and double-layer 
% potentials, in the case of circular scatterers.
% RCS = 10.*log10(2.*pi*(abs(Far_Field)).^2);
% -------------------------------------------------------------------------
% R = RCS(O, a, M_modes, k, theta, Density, Weight)
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat] : Vector of the centers of the scatterers
% a       [1 x N_scat] : Vector of the radius of the scatterers
% M_modes [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                        number of modes kept for the Fourier basis of the
%                        scatterer number "p".
% k            [1 x 1] : Wavenumber in the vacuum
% theta [1 x N_angles] : Angles of observation (rad).
%                        Note that theta could be a vector.
% Density      [N x 1] : Coefficients in the Fourier bases of the
%           or [N x 1]   density (or densities)
%                        N = 2*sum(M_modes+1).
% Weight       [1 x 2] : Weight to apply to the potentials (coef. alpha_p)
%      or [N_scat x 2]   if Weight is a line, then Weight(p,j) = Weight(j).
%
% Output Arguments :
% -----------------
% R     [N_angles x 1] : Vector of the Far Field F
%
%
%

function R = RCS(O, a, M_modes, k, theta, DensityCoef, TypeOfOperator)
    %Compute far field
    F = FarField(O, a, M_modes, k, theta, DensityCoef, TypeOfOperator);
    %Compute Radar Cross Section
    R = FarField_to_RCS(F);
end

