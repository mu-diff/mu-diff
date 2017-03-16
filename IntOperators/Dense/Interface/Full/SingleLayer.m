% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the matrix of the trace of the single layer integral operator 
% in the Fourier basis.
% -------------------------------------------------------------------------
% Trace of the single-layer operator L for a density rho in H^{-1/2}(Gamma)
% for all x in \Gamma :
%   L rho (x) = \int_{Gamma} G(x,y) rho(y) dy
%
% n = unit outward normal of the boundary Gamma of Omega, pointing
%     outside Omega
% G(x,y)=i/4*H_0^1(k\|\xx-\yy\|), zeroth order Hankel function of 1st kind
% -------------------------------------------------------------------------
%   A = SingleLayer(O, a, M_modes, k)
%
% Output Arguments :
% -----------------
% A [2*sum(M_modes+1) x 2*sum(M_modes+1)] : Matrix of the system
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat]  : Coordinates of the center of the disk 1
% a       [1 x N_scat]  : Radius of disk 1
% M_Modes [1 x N_scat]  : Truncation index in the Fourier series of
%                         obstacles
% k       [1 x 1]       : Wavenumber in the vacuum
% 
% See also BlockIntegralOperator, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator

%%
function A = SingleLayer(O, a, M_modes, k)
    A = IntegralOperator(O, a, M_modes, k, 2);
end
