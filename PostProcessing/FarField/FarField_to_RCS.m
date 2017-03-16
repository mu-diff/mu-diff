% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Given a Far field pattern, this function computes the Radar Cross
% Section (RCS).
% 
% Syntax
% ------
% R = FarField_to_RCS(Far_Field)
% 
% Input Arguments :
% -----------------
% Far_Field [N x 1]     :Far Field pattern.
%                        it is supposed to be a vector containing the
%                        far field pattern for every angle of [0,2*pi].
%                        N is the number of points used to discretize [0,2pi]
% Output Argument :
% -----------------
% rcs [N x 1]           :Radar Cross Section
%
% Result
% -------
% R = 10.*log10(2.*pi*(abs(Far_Field)).^2);
%

function R = FarField_to_RCS(Far_Field)
    R = 10.*log10(2.*pi*(abs(Far_Field)).^2);
end