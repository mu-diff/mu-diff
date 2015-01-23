% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Parser that change a character array into a scalar array
% 
% 
% 

function B = IntegralOperatorParser(A)
    ni = size(A,1);
    nj = size(A,2);
    nk = size(A,3);
    
    B = zeros(size(A));
    
    for i=1:ni
        for j =1:nj
            for j = 1:nk
                aijk = A(i,j,k);
                if(ischar(aijk))
                    switch aijk
                        case 'Z', bijk = 0;
                        case 'I', bijk = 1;
                        case 'L', bijk = 2;
                        case 'M', bijk = 3;
                        case 'N', bijk = 4;
                        case 'D', bijk = 5;
                        case 'Lsgl', bijk = 6;
                        case 'Dsgl', bijk = 7;
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