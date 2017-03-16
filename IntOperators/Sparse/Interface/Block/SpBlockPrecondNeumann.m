% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix of the integral operator obtained by
% preconditioning the dn double-layer by the single scattering of
% sound-hard obstacles (Neumann), from obstacle q to obstacle p, in the 
% Fourier basis
% SPARSE VERSION
% -------------------------------------------------------------------------
%   [LeftPart, MiddlePart, RightPart] = SpBlockPrecondNeumann(Op, ap, Np, Oq, aq, Nq, Nmax, k)
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

function [LeftPart, MiddlePart, RightPart] = SpBlockPrecondNeumann(Op, ap, Np, Oq, aq, Nq, Nmax, k)
    [LeftPart, MiddlePart, RightPart] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, k, 5);
end
