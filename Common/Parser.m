% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Parser that change a character array into a scalar array
% B = Parser(A)
%
% OUTPUT ARGUMENTS:
% -----------------
% B   [same size than A]       : Same matrix of A except that the eventual
%                                charactere has been replaced by integer
%
% INPUT ARGUMENTS:
% ----------------
% A      [? x ? x ?]    : Coordinates of the center of the disk p
%
% Table of correspondances:
% -------------------------
% Letter    Number      Operator
%'Z'        0           Null Matrix
%'I'        1           Identity
%'L'        2           Single Layer
%'M'        3           Double Layer
%'N'        4           Dn Single Layer
%'D'        5           Dn Double Layer
%'P'        6           Precond Dirichlet (single scat. precond)
%'Q'        7           Precond Neumann (single scat. precond)
% 
% See also IntegralOperator, BlockIntegralOperator

function B = Parser(A)
    ni = size(A,1);
    nj = size(A,2);
    nk = size(A,3);
    
    B = zeros(size(A));
    
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
                        case 'P', bijk = 6; % Precond Dirichlet
                        case 'Q', bijk = 7; % Precond Neumann
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