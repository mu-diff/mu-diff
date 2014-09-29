% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Returns A + alpha*I where A is a sparse operator already computed.
%   A = SpAddIdentity(A, alpha, M_modes)
% The addition with identity can be done directly in the creation of the
% matrix A. This function is hence used as a "post-processing". The
% alpha-identity is added to the "Left-Part" of A (ie: A{1}).
%
% See also SpIntegralOperator, SpBlockIntegralOperator, SpMatVec

function A = SpAddIdentity(A, alpha, M_modes)
    N_scat = size(A{1},2);
    for p = 1:N_scat
        Np = M_modes(p);
        A{1}(1:2*Np+1,p,p) = alpha + A{1}(1:2*Np+1,p,p);
    end
end