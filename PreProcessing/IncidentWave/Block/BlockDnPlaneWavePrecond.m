% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the vector of (the opposite of) the coefficients of the normal
% derivative trace of an incident plane wave on one obstacle for the 
% multiple scattering problem by disks, in the Fourier bases.
% The normal derivative trace itself is multiplied by the inverse of the 
% dn double-layer operator (single scattering preconditioner for sound hard
% obstacles)
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -u^inc
% -------------------------------------------------------------------------
%   Bp = BlockDnPlaneWavePrecond(Op, ap, Np, k, OS)
%
% OUTPUT ARGUMETNS:
% -----------------
% Bp [2*Np+1,1]   : Vector of the coefficients
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% Op          [2 x N_scat] : Vector of the center of the scatterer
% ap          [1 x N_scat] : Radius of the scatterer
% Np          [1 x N_scat] : Truncation index in the Fourier series of
%                            obstacles
% k           [1 x 1]      : Wavenumber in the vacuum
% beta_inc    [1 x 1]      : Direction of the wave
%
% See also IncidentWave, BlockIncidentWave, DnPlaneWavePrecond
%
%

function Bp = BlockDnPlaneWavePrecond(Op, ap, Np, k, beta_inc)
 Bp = BlockIncidentWave(Op, ap, Np, k, 6, beta_inc);
end