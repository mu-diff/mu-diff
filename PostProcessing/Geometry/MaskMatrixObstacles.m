% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Function that check whether a point is or not in a scatterer.
% This function return a matrix "res" of the same size of X or Y, obtained
% by "meshgrid".
% Coefficient res(i,j) is equal to:
%  p if the point (X(i,j), Y(i,j)) belongs to Omega_p
%  p+0.5 if the point is on the boundary of Omega_p
%  0 otherwise (outside every obstacles)
% The obstacles are supposed to be circular.
% ----------------------------------------------------------
% res = MaskMatrixObstacles(X, Y, O, a)
% 
% Input (N_scat = number of obstacle = length(a))
% -----
%   X [Nx x Ny]     : Meshgrid in X
%   Y [Nx x Ny]     : Meshgrid in Y
%   O [2 x N_scat]  : Matrix of the coordinates of the center of the
%                     scatterers
%   a [1 x N_scat]  : Vector of the radii of the scatterers
%
% Output
% ------
%   res [Nx x Ny]   : Matrix where the coefficient (i,j) is equal to:
%                     0 if the point X(i,j), Y(i,j) is outside the scatterers
%                     p if the point X(i,j), Y(i,j) is inside the p^th scatterer
%                     p+0.5 if the point X(i,j), Y(i,j) is on the boundary 
%                     of the p^h scatterers
%
%%
function res = MaskMatrixObstacles(X, Y, O, a)

tolM = 0;

N_scat = length(a);
res = zeros(size(X));

for p=1:N_scat
   ap = a(p);
   Op = O(:,p);
   rp = sqrt((X-Op(1)).^2 + (Y-Op(2)).^2);
   res(rp < ap-tolM) = p; % setting "p" to the points belonging to obstacle p
   
   %boundary problem
   rp_minus = find((rp >= ap-tolM));
   rp_plus = find((rp <= ap+tolM));   
   res(intersect(rp_plus, rp_minus)) = p + 0.5; % setting "p+0.5" on the boundary
end