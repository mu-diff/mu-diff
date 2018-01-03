% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Sparse matrix vector product between possibly multiple matrices A and a 
% vector X. Matrices are stored as a cell and have the special structure:
% A{1} = LeftPart  (off-diagonal blocks + diagonal blocks)
% A{2} = Toeplitz Part (off-diagonal blocks)
% A{3} = RightPart  (off-diagonal blocks)
% -------------------------------------------------------------
% The matrix vector product is done through a cross-correlation using xcorr
% Matlab function.
% -------------------------------------------------------------
%   Y = SpMatVec(X, M_modes, ListOfOperators, varargin)
%
% OUTPUT ARGUMENT:
% ----------------
% Y [size(X)] : A*X
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------------------------------------
%
% X       [sum(2*M_modes+1), 1] : Vector
% M_Modes [N_scat, 1]           : Truncation index in the Fourier series of
%                                 obstacles
% ListOfOperators  cell(NOp,1)  : Cell array of integral operators (sparse storage)
%
% OPTION (weight):
% ----------------
% Y = SpMatVec(X, M_modes, ListOfOperators, ListOfWeight)
%
% ListOfWeight [length(ListOfOperators), 1] : list of multiplicative constant 
%                                             to apply to the operators
%
% EXAMPLE:
% --------
% Compute Y =(0.5*I + N)X, where I = identity and N = DnSingleLayer:
%   > [A{1},A{2},A{3}] = SpDnSingleLayer(O, a, M_modes, k)
%   > [I{1},I{2},I{3}] = SpIdentity(O, a, M_modes);
%   > Y = SpMatVec(X, M_modes, [I,A], [0.5, 1]);
%
% See also SpSingleMatVec, SpBlockIntegralOperator, SpIntegralOperator, xcorr
%


function Y = SpMatVec(X, M_modes, ListOfOperators, varargin)

noperator = size(ListOfOperators, 2);

if(length(varargin) >= 1)
    Weight = varargin{1};
    if(iscell(Weight))
        Weight = cell2mat(Weight);
    end
else
    Weight = ones(1,noperator);
end

%To do
if(length(varargin) >= 2)
    OpMask = varargin{2};
else
    OpMask = ones(1,noperator-1);
end

%Check
if(length(Weight) ~= noperator)
   error('Weight must be of size size(ListOfOperators,2)'); 
end
if(length(OpMask) ~= noperator-1)
   error('There must be ''size(ListOfOperators,2)-1'' operations'); 
end

Y = zeros(size(X));
for iOper = 1:noperator
    if(noperator == 1)
        A = ListOfOperators;
    else
        A = ListOfOperators{iOper};
    end
    this_weight = Weight(iOper);
    Y = Y + this_weight*SpSingleMatVec(X, M_modes, A);
    clear A this_weight;
end
