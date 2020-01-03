% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix of the normal derivative trace of the 
% single layer integral operator from obstacle q to obstacle p, in the 
% Fourier basis:
%       (N_q rho_q) |_{\Gamma_p}
% IMPORTANT REMARK: the jump relation are NOT taken into account and the 
% identity operator is thus NOT contained in the obtained matrix. So it is 
% not "really" the normal trace of the single layer if Gammap == Gammaq.
% -------------------------------------------------------------------------
% Normal derivative trace of the Single-Layer operator for a density 
% rho_q in H^{-1/2}(Gamma_q)
% for all x in \Gamma_p :
%   N_q rho_q (x) = \int_{Gamma_q} dn_x G(x,y) rho_q(y) dy
%
% n = unit outward normal of the boundary Gamma_p of Omega_p, pointing
%     outside Omega_p
% G(x,y)=i/4*H_0^1(k\|\xx-\yy\|), zeroth order Hankel function of 1st kind
% -------------------------------------------------------------------------
%   Apq = BlockDnSingleLayer(Op, ap, Np, Oq, aq, Np, k)
%
% Output Arguments :
% -----------------
% Apq [2*Np+1 x 2*Nq+2] : Matrix of the system
%
% Input Arguments :
% -----------------
% Op   [2 x 1]  : Coordinates of the center of the disk p
% ap   [1 x 1]  : Radius of disk p
% Np   [1 x 1]  : Truncation index in the Fourier series of obstacle p
% Oq   [2 x 1]  : Coordinates of the center of the disk q
% aq   [1 x 1]  : Radius of disk q
% Nq   [1 x 1]  : Truncation index in the Fourier series of obstacle q
% k    [1 x 1]  : Wavenumber in the vacuum
% 
% See also BlockIntegralOperator, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator

%%
function Apq = BlockDnSingleLayer(Op, ap, Np, Oq, aq, Nq, k)
    Apq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 4);
end
