% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the derivatives of the Neumann function Y
% (Or Bessel function of second kind)
% dy=dbessely(m,x) 
% m = order, x = evaluation point  
% Taken from the help of bessely:
%    If m and x are arrays of the same size, the result is also that size.
%    If either input is a scalar, it is expanded to the other input's size.
%    If one input is a row vector and the other is a column vector, the
%    result is a two-dimensional table of function values.
function dy = dbessely(m,k)

dy = (bessely(m-1,k)-bessely(m+1,k))/2;