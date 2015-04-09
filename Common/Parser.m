% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Parser that change a character or numeric (or mix) array or cell 
% into a scalar array of corresponding integral operators.
% B = Parser(A)
%
% OUTPUT ARGUMENTS:
% -----------------
% B   [same size than A]       : Array or numerical values only.
%
% INPUT ARGUMENTS:
% ----------------
% A : can be a one, two or three dimensionnal array or cell.
% If A is an array, it can be composed by numeric or (only one) char values
% If A is a cell, it can be composed by numeric or char values (more than 1
% char). See below.
%
% Table of correspondances:
% -------------------------
% Letter (cell case)   Number    Associated Operator
%'Z'                     0       Null Matrix
%'I'                     1       Identity
%'L'                     2       Single Layer
%'M'                     3       Double Layer
%'N'                     4       Dn Single Layer
%'D'                     5       Dn Double Layer
%'P' ('Lprec')           6       Precond Dirichlet (single scat. precond, diag(L)^-1*L)
%'Q' ('Dprec')           7       Precond Neumann (single scat. precond, diag(D)^-1*D)
% 
% See also IntegralOperator, BlockIntegralOperator, SpIntegralOperator, 
% SpBlockIntegralOperator

function B = Parser(A)
    ni = size(A,1);
    nj = size(A,2);
    nk = size(A,3);
    
    B = zeros(size(A));

    if(iscell(A))
        for i=1:ni
            for j =1:nj
                for k = 1:nk
                    aijk = A{i,j,k};
                    if(ischar(aijk))
                        switch aijk
                            case 'Z', bijk = 0; % Null
                            case 'I', bijk = 1; % Identity
                            case 'L', bijk = 2; % Single Layer
                            case 'M', bijk = 3; % Double Layer
                            case 'N', bijk = 4; % Dn Single Layer
                            case 'D', bijk = 5; % Dn Double Layer
                            case 'Lprec', bijk = 6; % Precond Dirichlet (diag(L)^-1*L)
                            case 'Dprec', bijk = 7; % Precond Neumann (diag(D)^-1*D)
                            case 'P', bijk = 6; % Precond Dirichlet (diag(L)^-1*L)
                            case 'Q', bijk = 7; % Precond Neumann (diag(D)^-1*D)
                            otherwise, error(['Parser error: unknown value ', aijk]);
                        end
                    elseif(isscalar(aijk))
                        bijk = aijk;
                    else
                        error('Parser error: unknown type (non scalar nor char)');
                    end
                    B(i,j,k) = bijk;
                end
            end
        end
    else %mixte of numeric and char
        for i=1:ni
            for j =1:nj
                for k = 1:nk
                    aijk = A(i,j,k);
                    if(ischar(aijk))
                        switch aijk
                            case 'Z', bijk = 0; % Null
                            case 'I', bijk = 1; % Identity
                            case 'L', bijk = 2; % Single Layer
                            case 'M', bijk = 3; % Double Layer
                            case 'N', bijk = 4; % Dn Single Layer
                            case 'D', bijk = 5; % Dn Double Layer
                            case 'P', bijk = 6; % Precond Dirichlet (diag(L)^-1*L)
                            case 'Q', bijk = 7; % Precond Neumann (diag(D)^-1*D)
                            otherwise, error(['Parser error: unknown value ', aijk]);
                        end
                    elseif(isscalar(aijk))
                        bijk = aijk;
                    else
                        error('Parser error: unknown type (non scalar nor char)');
                    end
                    B(i,j,k) = bijk;
                end
            end
        end
    end
end