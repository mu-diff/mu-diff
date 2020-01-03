% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the vector of (the opposite of) the coefficients of the trace
% of the Green function for the multiple scattering problem by disks, in 
% the Fourier bases.
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -G|_\Gamma
% -------------------------------------------------------------------------
% Green function: G(x,y) = 0.25*i*H_0^1(k\|x-y\|)
% -------------------------------------------------------------------------
%   B = PointSource(O, a, M_modes, k, OS)
%
% OUTPUT ARGUMETNS:
% -----------------
% B [sum(2*M_modes+1),1]    : Vector of the coefficients
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% O           [2 x N_scat] : Vector of the centers of the scatterers
% a           [1 x N_scat] : Vector of the radius of the scatterers
% M_Modes     [1 x N_scat] : Truncation index in the Fourier series of
%                            obstacles
% k           [1 x 1]      : Wavenumber in the vacuum
% OS          [2 x 1]      : Coordinates of the point source
%
% See also IncidentWave, PlaneWave, DnPlaneWave, DnPointSource,
% PlaneWavePrecond, DnPlaneWavePrecond
%
%

function B = PointSource(O, a, M_modes, k, OS)
    B = IncidentWave(O, a, M_modes, k, 3, OS);
end