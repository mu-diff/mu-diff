% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Parser that change a character or numeric (or mix) array or cell 
% into a scalar array of corresponding incident waves.
% B = ParserIncidentWave(A)
%
% OUTPUT ARGUMENTS:
% -----------------
% B   [same size than A]       : Array or numerical values only.
%
% INPUT ARGUMENTS:
% ----------------
% A : can be a one or two dimensionnal array or cell.
% If A is an array, it can be composed by numeric values only
% If A is a cell, it can be composed by numeric or char values (more than 1
% char). See below.
%
% Table of correspondances:
% -------------------------
%    Name                   Number
% 'PlaneWave'                 1
% 'DnPlaneWave'               2
% 'PointSource'               3
% 'DnPointSource'             4
% 'PlaneWavePrecond'          5
% 'DnPlaneWavePrecond'        6
% 'PointSourcePrecond'        7
% 'DnPointSourcePrecond'      8
% 
% See also IncidentWave, BlockIncidentWave

function B = ParserIncidentWave(A)
    if (ischar(A))
        B = ParseCharIncidentWave(A);
    else
        nrows = size(A,1);
        B = zeros(size(A));
        if(iscell(A))
            for i=1:nrows
                if(iscell(A))
                    ai = A{i};
                else
                   ai = A;
                end
                if(ischar(ai))
                    bi = ParseCharIncidentWave(ai);
                elseif(isscalar(ai))
                    bi = ai;
                else
                    error('ParserIncidentWave error: unknown type (non scalar nor char)');
                end
                B(i) = bi;
            end
        else %one value only
            B = A;
        end
    end
end


%%
function bi = ParseCharIncidentWave(ai)
    switch ai
        case 'PlaneWave', bi = 1;
        case 'DnPlaneWave', bi = 2;
        case 'PointSource', bi = 3;
        case 'DnPointSource', bi = 4;
        case 'PlaneWavePrecond', bi = 5;
        case 'DnPlaneWavePrecond', bi = 6;
        case 'PointSourcePrecond', bi = 7;
        case 'DnPointSourcePrecond', bi = 8;
        otherwise, error(['ParserIncidentWave error: unknown value ', ai]);
    end
end