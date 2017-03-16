% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% function which being given two points O1 and O2 computes the angle
%% between O2O1 and the unit vector along the abscissa OI

function alpha = fangle(O1,O2)
m_1 = O1(1); n_1 = O1(2);
m_2 = O2(1); n_2 = O2(2);
O2O1= [m_1-m_2, n_1-n_2];
OI=[1,0];

d12 = norm(O2O1);
if d12 == 0
    alpha = 0;
else  

if    det([OI;O2O1]) >= 0
    alpha = acos((m_1-m_2)/d12);
   
 else 
    alpha=-acos((m_1-m_2)/d12) +2*pi;
    
end
end