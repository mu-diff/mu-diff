% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix Cpq of the Calderon projector C given by
% C = [I/2+M, L; D, I/2+N]
% and hence
% Cpq = [0.5*Ipq+Mpq, Lpq ; Dpq, 0.5*Ipq+Npq]
% Note that Ipq = 0 if p~=q and Ipq = identity otherwise
% -------------------------------------------------------------------------
%   Apq = BlockCalderonProjector(Op, ap, Np, Oq, aq, Nq, k)
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
% See also CalderonProjector, BlockIntegralOperator, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator

%%
function Apq = BlockCalderonProjector(Op, ap, Np, Oq, aq, Nq, k)
    Ipq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 1);
    Lpq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 2);
    Mpq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 3);
    Npq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 4);
    Dpq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, 5);
    Apq = [0.5*Ipq+Mpq, Lpq; Dpq, 0.5*Ipq+Npq];
end
