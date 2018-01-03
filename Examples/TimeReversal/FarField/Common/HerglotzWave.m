% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Computes the Herglotz wave with coefficient V and on the angles of emission
% from EmissionAngles. The computation is made on a grid (X,Y) with a
% wavenumber k. The integral is approximated using step functions with
% step halpha.

function U = HerglotzWave(V, k, halpha, EmissionAngles, X,Y)

U = zeros(size(X));

cosalpha = cos(EmissionAngles);
sinalpha = sin(EmissionAngles);
Nb_angles = length(EmissionAngles);

for n = 1: Nb_angles
    U = U + halpha * V(n)*exp(1i*k*(cosalpha(n)*X + sinalpha(n)*Y));
end
