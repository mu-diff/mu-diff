% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the vector of (the opposite of) the coefficients of the trace of
% an incident plane wave for the multiple scattering problem by disks, in 
% the Fourier bases.
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -u^inc
% -------------------------------------------------------------------------
% Trace of an incident plane wave of direction Beta (angle in radian)
% on the boundary Gamma :
% u^{inc}(x) = exp(ik(cos(Beta).x_1 + sin(Beta).x_2)),    for all x=(x_1,x_2) in Gamma
% -------------------------------------------------------------------------
%   B = IncidentWave(O, a, M_modes, k, beta_inc)
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
% See also IncidentWave, DnPlaneWave, PointSource, DnPointSource,
% PlaneWavePrecond, DnPlaneWavePrecond
%
%

function B = PlaneWave(O, a, M_modes, k, beta_inc)
 B = IncidentWave(O, a, M_modes, k, 1, beta_inc);
end