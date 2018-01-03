% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
function DiagCell = CellDiag(StringValue, varargin)
    DiagCell = cell(varargin{1:end});
    ni = size(DiagCell,1);
    nj = size(DiagCell,2);
    nk = size(DiagCell,3);
    if( nk > 1)
        error('cellDiag is for 2d use only');
    end
    mindim = min(ni,nj);
    for i = 1:mindim
       DiagCell{i,i} = StringValue; 
    end
end