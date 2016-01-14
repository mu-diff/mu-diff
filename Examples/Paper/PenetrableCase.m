% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
%
% Contributed by X. Antoine and B. Thierry
%
% Benchmark file of solving multiple scattering problem of an acoustic
% incident wave by a collection of homogeneous penetrable circular obstacles.
% --------------------------------------------------------------------
% The Radar Cross Section (RCS) is computed together with a map of the
% total field
% --------------------------------------------------------------------
% The integral formulation is based on a single-layer potential only
% --------------------------------------------------------------------
% Obstacles have an electric permittivity and magnetic permeability
% given by (resp.) eps_m and mu_m
% Ouside, elect. perm. is equal to eps_0 and magn. perm. to mu_0
% (see below)

disp('-----------------------');
disp('    New simulation     ');
disp('-----------------------');

clear all;
close all;
%% -----------------------------------------
%%            Params
%% -----------------------------------------
ind_fig = 1;
%% colors
ColorBIE = 'k-';
ColorBIE_Precond = 'k-.';

%% -----------------------------------------
%%            Pre Processing
%% -----------------------------------------
%% Incident Wave (point source)
XS = [0; 0];
%% GRID on which the fiels are computed (for the end)
XXmin = -30; XXmax = 30;
YYmin = -30; YYmax = 30;
lc = min(2*pi/15, 1/10); % characteristic length

%% Geometry (chose one)
bx = 2.1; by = 2.1;
Nx = 21; Ny = 21;
O_old = RectangularLattice(bx, by, Nx, Ny, 'Center', [0,0]);
a_old = 1*ones(1,size(O_old,2));
%Remove the line of disk with either X==0 or Y==0
[O,a] = RemoveDisk(O_old, a_old, 'X', 0, 'Y', 0);
N_scat = size(O,2);
%% Physical properties
% Wavenumbers
k = 1;
k_int = 2*k;
contrast = k_int./k;

%% Fourier series: truncation indices
tolM = 10.^(-10);
M_modes_plus = FourierTruncation(a, k, 'Min', 10);
M_modes_minus = FourierTruncation(a, k_int, 'Min', 10);
% Truncation indices are chosen as the max between intern and extern values
M_modes = max(M_modes_plus, M_modes_minus);
sum_modes = sum(2*M_modes+1);
%% GMRES parameters
RESTART = [];
TOL = 10^(-10);
MAXIT = 500;

%% -----------------------------------------
%            Integral formulation
% Single-Layer potential (plus and minus sign refer to ouside and
% inside the obstacles)
% u = \int_{\Gamma} G^{+/-}(x,y) rho^{+/-}(y) d\Gamma(y),  x
% rho^{+/-} = densities (unknown)
% vector of the "lambda" (contrast)

%% Boundary integral operators
%Single-layer (exterior)
Splus = SingleLayer(O, a, M_modes, k); 
%Single-layer (intern: no off-diagonal blocks but different wave numbers)
Sminus = IntegralOperator(O, a, M_modes, k_int, 2*eye(N_scat,N_scat));
%Dn Single Layer (exterior)
Nplus = DnSingleLayer(O, a, M_modes, k);
%Dn Single Layer (intern, no off-diagonal blocks but different wave numbers)
Nminus = IntegralOperator(O, a, M_modes, k_int, 4*eye(N_scat,N_scat));
%identity matrix
Identity = eye(size(Nplus));
%% Boundary integral operator
A = [Splus, -Sminus; -0.5*Identity + Nplus, -(0.5*Identity + Nminus)];
%cleaning memory
clear Splus Sminus Nplus Nminus Identity ;
%% Incident plane wave (right hand side)
Uinc = PointSource(O, a, M_modes, k, XS);
DnUinc = DnPointSource(O, a, M_modes, k, XS);
[B] = [Uinc;DnUinc];
clear Uinc DnUinc;
%% Solving linear system
rho = A\B;
%% Extracting the "exterior" and "interior" densities
rho_plus = rho(1:sum_modes);
rho_minus = rho(sum_modes+1:end);
clear rho;
%% Radar Cross Section (RCS)
theta_RCS = 0:360;
theta_RCS_rad = theta_RCS*2*pi/360;

R = RCS(O, a, M_modes, k, theta_RCS_rad, rho_plus, [1,0]);
%% Scattered field on grid
%SEE BEGINING FOR XXmin, XXmax, YYmin and YYmax values
XX = [XXmin:lc:XXmax];
YY = [YYmin:lc:YYmax];
[X,Y] = meshgrid(XX,YY);

Ue = ExternalPotential(X, Y, O, a, M_modes, k, rho_plus, [1,0]);
Ui = InternalPotential(X, Y, O, a, M_modes, k_int, rho_minus, [1,0], 'OnBoundary', 1);
U = Ue + Ui;

%% Incident wave on grid (outside obstacles)
UincOnMesh = IncidentWaveOnGrid(X, Y, k, 'PointSource', XS);
Matrix_Not_Obstacles = (MaskMatrixObstacles(X, Y, O, a) == 0);
UincOnMesh = UincOnMesh.*Matrix_Not_Obstacles;

%% Total field
U_tot = U + UincOnMesh;
%% Save
clear A B;
save('data.mat');