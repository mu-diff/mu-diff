% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the total, scattered and incident waves on a grid, for the 
% penetrable problem of multiple scattering by disk. The solution must have 
% been computed with PenetrableSolver.
% -------------------------------------------------------------------------
%   [U_tot, U, UincOnMesh] = PenetrableNearField(Solution, X, Y)
%
% OUTPUT ARGUMENTS:
% -----------------
% U_tot       (sizeof(X)) : Total field
% U           (sizeof(X)) : Scattered field
% UincOnMesh  (sizeof(X)) : Incident wave
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by PenetrableSolver
% X         [Nx x Ny]   : X-Grid point
% Y         [Nx x Ny]   : Y-Grid point
%
% See also PenetrableSolver, PenetrableFarField, PenetrableRCS, meshgrid
%
function [U_tot, U, UincOnMesh] = PenetrableNearField(Solution, X, Y)

O = Solution{1};
a = Solution{2};
k = Solution{3};
k_int = Solution{4};
M_modes = Solution{5};
TypeOfWave = Solution{6};
ParamWave = Solution{7};
rho_plus = Solution{8};
rho_minus = Solution{9};

%% Scattered field
Ue = ExternalPotential(X, Y, O, a, M_modes, k, rho_plus, [1,0]);
Ui = InternalPotential(X, Y, O, a, M_modes, k_int, rho_minus, [1,0], 'OnBoundary', 1);
U = Ue + Ui;
clear Ue Ui;
%% Incident wave
UincOnMesh = IncidentWaveOnGrid(X, Y, k, TypeOfWave, ParamWave);
Matrix_Not_Obstacles = MaskMatrixObstacles(X, Y, O, a) == 0;
UincOnMesh = UincOnMesh.*Matrix_Not_Obstacles; %delete points in obstacles
%% Total field
U_tot = U + UincOnMesh;
