% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the identity matrix.
% -------------------------------------------------------------------------
%   I = Identity(O, a, M_modes)
%
% Output Arguments :
% -----------------
% I [2*sum(M_modes+1) x 2*sum(M_modes+1)] : Matrix of the system
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat]  : Coordinates of the center of the disk 1
% a       [1 x N_scat]  : Radius of disk 1
% M_Modes [1 x N_scat]  : Truncation index in the Fourier series of
%                         obstacles
% 
% See also BlockIntegralOperator, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator, SpIdentity

%%
function I = Identity(O, a, M_modes)
    I = IntegralOperator(O, a, M_modes, 1, 1);
end
