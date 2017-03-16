% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% From an array MaskArray, build a cell of the exact same size and
% plug the string (or cell of strings) StringValue at the position where
% MaskMatrix has non zero coefficient.
% -------------------------------------------------------------------------
% function ouputCell = CellFromMat(StringValue, MaskArray)
% -------------------------------------------------------------------------
% OUTPUT ARGUMENTS:
% -----------------
% ouputCell  {Nx x Ny x ...}  : ouput cell array
%
% INPUT ARGUMENTS:
% ----------------
% StringValue     char     : Char to be placed in output cell
%         or cell of char
% MaskArray [Nx x Ny x ...]: Mask array
%
% See also cell

function ouputCell = CellFromMat(StringValue, MaskArray)
    ouputCell = cell(size(MaskArray));
    NonZeroList = find(MaskArray ~= 0);
    nNonZero = length(NonZeroList);
    for i=1:nNonZero
        ouputCell{NonZeroList(i)} = StringValue;
    end
end