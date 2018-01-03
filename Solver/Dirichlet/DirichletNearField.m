% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the total, scattered and incident waves on a grid, for the 
% Dirichlet problem of multiple scattering by disk. The solution must have 
% been computed with DirichletSolver.
% -------------------------------------------------------------------------
%   [U_tot, U, UincOnMesh] = DirichletNearField(Solution, X, Y)
%
% OUTPUT ARGUMENTS:
% -----------------
% U_tot       (sizeof(X)) : Total field
% U           (sizeof(X)) : Scattered field
% UincOnMesh  (sizeof(X)) : Incident wave
%
% INPUT ARGUMENTS
% ---------------
% Solution   cell(7,1)  : Solution computed by DirichletSolver
% X         [Nx x Ny]   : X-Grid point
% Y         [Nx x Ny]   : Y-Grid point
%
% See also DirichletSolver, DirichletFarField, DirichletRCS, meshgrid
%
function [U_tot, U, UincOnMesh] = DirichletNearField(Solution, X, Y)

O = Solution{1};
a = Solution{2};
k = Solution{3};
M_modes = Solution{4};
TypeOfWave = Solution{5};
ParamWave = Solution{6};
rho = Solution{7};

%% Scattered field
U = ExternalPotential(X, Y, O, a, M_modes, k, rho, [1,0], 'OnBoundary', 1);
%% Incident wave
UincOnMesh = IncidentWaveOnGrid(X, Y, k, TypeOfWave, ParamWave);
Matrix_Not_Obstacles = MaskMatrixObstacles(X, Y, O, a) == 0;
UincOnMesh = UincOnMesh.*Matrix_Not_Obstacles; %delete points in obstacles
%% Total field
U_tot = U + UincOnMesh;
