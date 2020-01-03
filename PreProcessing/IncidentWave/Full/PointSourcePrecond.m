% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the vector of (the opposite of) the coefficients of the trace of
% an incident point source wave for the multiple scattering problem by disks,
% in the Fourier bases.
% The trace itself is multiplied by the inverse of the single-layer
% operator (single scattering preconditioner for sound soft obstacles)
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -u^inc
% -------------------------------------------------------------------------
%   B = PointSourcePrecond(O, a, M_modes, k, beta_inc)
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
% beta_inc    [1 x 1]      : Direction of the wave
%
% See also IncidentWave, PlaneWave, DnPlaneWave, PointSource, DnPointSource,
% PlaneWavePrecond, DnPlaneWavePrecond, DnPointSourcePrecond
%
%

function B = PointSourcePrecond(O, a, M_modes, k, beta_inc)
 B = IncidentWave(O, a, M_modes, k, 7, beta_inc);
end