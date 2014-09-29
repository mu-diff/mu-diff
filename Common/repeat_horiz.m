% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Repeat horizontally N_repeat times a matrix or a vector
% 
% Syntax
% ------
% Res = repeat_horiz(mat, N_repeat)
%
% Input Arguments
% ---------------
% mat      [M x N]      : Matrix to be repeated
% N_repeat [1x1]        : Number of repetitions
% 
% Example
% -------
% mat = [1,2,3;4,5,6];
% res = repeat_horiz(mat,3);
% res =
%     1     2     3     1     2     3     1     2     3
%     4     5     6     4     5     6     4     5     6

function res = repeat_horiz(mat, N_repeat)

    res = [];
    for cpt = 1:N_repeat
       res = [res,mat];
    end
end