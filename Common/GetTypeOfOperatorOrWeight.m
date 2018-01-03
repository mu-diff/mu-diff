% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
%
%

%%
function value = GetTypeOfOperatorOrWeight(A, p, q)
    if (isrow(A) || iscolumn(A))
       value = A;
    elseif(size(A,3) == 1)
       value = A(p,q);
    else
       value = A(:,p,q);
    end
end
