% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the derivatives of the Hankel function H_m^{(K)} of the first
% kind (K =1) or the second kind( K =2).
% dh = dbesselh(m,K,x)
% m = order, K = kind (1 or 2), x = evaluation point  
% Taken from the help of bessely:
%    If m and x are arrays of the same size, the result is also that size.
%    If either input is a scalar, it is expanded to the other input's size.
%    If one input is a row vector and the other is a column vector, the
%    result is a two-dimensional table of function values.
% Note : this function involves dbesselj and dbessely.

function dh = dbesselh(m,K,x)

    if (K == 1)
     dh = dbesselj(m,x) + 1i.*dbessely(m,x);
    elseif( K == 2)
     dh = dbesselj(m,x) - 1i.*dbessely(m,x);
    else
        error('K is not equal to 1 or 2');
    end