% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix of the integral operator obtained by
% preconditioning the single-layer by the single scattering of sound-soft
% obstacles (Dirichlet), from obstacle q to obstacle p, in the Fourier basis
% -------------------------------------------------------------------------
%   Apq = BlockPrecondDirichlet(Op, ap, Np, Oq, aq, Np, k)
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
function Apq = BlockPrecondDirichlet(Op, ap, Np, Oq, aq, Nq, k)
    Apq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 6);
end
