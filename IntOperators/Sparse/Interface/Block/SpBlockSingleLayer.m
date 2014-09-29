% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix of the trace of the single layer integral 
% operator from obstacle q to obstacle p, in the Fourier basis:
%       (L_q rho_q) |_{\Gamma_p}
% SPARSE VERSION
% -------------------------------------------------------------------------
% Trace of the Single-Layer operator for a density rho_q in H^{-1/2}(Gamma_q)
% for all x in \Gamma_p :
%   L_q rho_q (x) = \int_{Gamma_q} G(x,y) rho_q(y) dy
%
% n = unit outward normal of the boundary Gamma_p of Omega_p, pointing
%     outside Omega_p
% G(x,y)=i/4*H_0^1(k\|\xx-\yy\|), zeroth order Hankel function of 1st kind
% -------------------------------------------------------------------------
%   [LeftPart, MiddlePart, RightPart] = SpBlockSingleLayer(Op, ap, Np, Oq, aq, Nq, Nmax, k)
%
% OUTPUT ARGUMENTS:
% -----------------
% LeftPart   [2*Nmax+1, 1]       : Left part in the off-diag. blocks and
%                                  also the diagonal blocks of the matrix
% MiddlePart [2*(2*Nmax+1)-1, 1] : Middle part (Toeplitz part) of the 
%                                  off-diag. blocks
% RightPart  [2*Nmax+1, 1]       : Right part of the off-diag. blocks
%
%
% INPUT ARGUMENTS:
% ----------------
% Op              [2 x 1]    : Coordinates of the center of the disk p
% ap              [1 x 1]    : Radius of disk p
% Np              [1 x 1]    : Truncation index in the Fourier series of
%                              obstacle p
% Oq              [2 x 1]    : Coordinates of the center of the disk q
% aq              [1 x 1]    : Radius of disk q
% Nq              [1 x 1]    : Truncation index in the Fourier series of 
%                              obstacle q
% Nmax            [1 x 1]    : Reference number N (max of every other Np)
% k               [1 x 1]    : Wavenumber in the vacuum
% 
% See also SpBlockIntegralOperator, SpIntegralOperator, SpMatVec,
% BlockIntegralOperator, IntegralOperator
% 

%%
function [LeftPart, MiddlePart, RightPart] = SpBlockSingleLayer(Op, ap, Np, Oq, aq, Nq, Nmax, k)
    [LeftPart, MiddlePart, RightPart] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, k, 0);
end
