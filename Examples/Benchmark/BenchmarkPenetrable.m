% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
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
%% Incident Wave (chose one)
%IncidentWave = 'PlaneWave';
IncidentWave = 'PointSource';
% Incidence angle
beta_inc = pi;
%For point source
XS = [-10; 0];
%% GRID on which the fiels are computed (for the end)
XXmin = -10; XXmax = 10;
YYmin = -10; YYmax = 10;
lc = min(2*pi/15, 1/10); % characteristic length

%% Geometry (chose one)
GEOMETRY = 'manual';
%GEOMETRY = 'random';
%GEOMETRY = 'rectangular';
%GEOMETRY = 'triangular';

switch GEOMETRY
    case 'manual',
        O = [-2.5, 2.5    ; -2, 1];
        a = [1.5, 0.75];
    case 'random',
        x_min = -3; x_max = 3;
        y_min = -3; y_max = 3;
        distance_min = 0.01;
        a_min = 0.09;
        a_max = 0.11;
        N_scat = 3;
        [O,a] = CreateRandomDisks(x_min,x_max,y_min,y_max, N_scat, a_min, a_max, distance_min);
    case 'rectangular',        
        bx = 3; by = 3;
        Nx = 2; Ny = 2;
        O = RectangularLattice(bx, by, Nx, Ny, 'Center', [-2,-2]);
        a = 1*ones(1,size(O,2));
    case 'triangular',
        bx = 1; by = 1;
        Nx = 2; Ny = 2;
        O = TriangularLattice(bx, by, Nx, Ny, 'Center', [-2,-2]);
        a = 0.1726*ones(1,size(O,2));
    otherwise,
        error('Please, chose a geometry');
end
N_scat = size(O,2);
%% Drawing the circular obstacles
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2);
%axis([x_min, x_max, y_min, y_max]);
ind_fig = ind_fig +1;
xlabel('x'); ylabel('y');
title('Obstacles');
axis equal;
%% Physical properties
omega    = 2*pi; %pulsation
%Vaccum:
epsilon0 = 1/omega;
mu0      = 1/omega;
%Obstacles:
epsilon_m = 2.2^2*epsilon0;
mu_m      = mu0;

% Checking for errors:
if(length(epsilon_m) < N_scat)
    epsilon_m = [epsilon_m, epsilon_m(end)*ones(1, N_scat-length(epsilon_m))];
elseif(length(epsilon_m) > N_scat)
    epsilon_m = epsilon_m(1:N_scat);
end
if(length(mu_m) < N_scat)
    mu_m = [mu_m, mu_m(end)*ones(1, N_scat-length(mu_m))];
elseif(length(mu_m) > N_scat)
    mu_m = mu_m(1:N_scat);
end

% Wavenumbers
k = omega*sqrt(epsilon0 * mu0);  %exterior wavenumber
k_int = omega .* sqrt(epsilon_m .* mu_m);
contrast = k_int./k;

%% Fourier series: truncation indices
tolM = 10.^(-10);
M_modes_plus = FourierTruncation(a, k, 'Min', 4);
M_modes_minus = FourierTruncation(a, k_int, 'Min', 4);
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
%% -----------------------------------------
mu_matrix = zeros(sum_modes,sum_modes);
Sp = 0;
for p = 1:N_scat
    Np = 2*M_modes(p)+1;
    mu_matrix(Sp+1:Sp+Np, Sp+1:Sp+Np) = (mu0/mu_m(p))*eye(Np);
    Sp = Sp +Np;
end

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
A = [Splus, -Sminus; -0.5*Identity + Nplus, -mu_matrix*(0.5*Identity + Nminus)];
%cleaning memory
clear Splus Sminus Nplus Nminus Identity mu_matrix;
%% Incident plane wave (right hand side)
switch IncidentWave
    case 'PlaneWave',
        Uinc = PlaneWave(O, a, M_modes, k, beta_inc);
        DnUinc = DnPlaneWave(O, a, M_modes, k, beta_inc);
    case 'PointSource',
        Uinc = PointSource(O, a, M_modes, k, XS);
        DnUinc = DnPointSource(O, a, M_modes, k, XS);
end
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

ind_fig = ind_fig +1; figure(ind_fig);
plot(theta_RCS, R, ColorBIE);
title('Radar Cross Section');
xlabel('Angle of reception (degree)');
ylabel('Radar Cross Section (dB)');
axis tight;

%% Scattered field on grid
%SEE BEGINING FOR XXmin, XXmax, YYmin and YYmax values
XX = [XXmin:lc:XXmax];
YY = [YYmin:lc:YYmax];
[X,Y] = meshgrid(XX,YY);

Ue = ExternalPotential(X, Y, O, a, M_modes, k, rho_plus, [1,0]);
Ui = InternalPotential(X, Y, O, a, M_modes, k_int, rho_minus, [1,0], 'OnBoundary', 1);
U = Ue + Ui;

%% Incident wave on grid (outside obstacles)
switch IncidentWave
    case 'PlaneWave',
        UincOnMesh = IncidentWaveOnGrid(X, Y, k, 'PlaneWave', beta_inc);
    case 'PointSource',
        UincOnMesh = IncidentWaveOnGrid(X, Y, k, 'PointSource', XS);
end
Matrix_Not_Obstacles = (MaskMatrixObstacles(X, Y, O, a) == 0);
UincOnMesh = UincOnMesh.*Matrix_Not_Obstacles;

%% Total field
U_tot = U + UincOnMesh;

%% Display on mesh
% Note: to erase white artefacts, use
% set(gcf,'Renderer','Zbuffer');

ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, real(U_tot));
shading interp;
title(['Real part of the total field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off
 
ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, imag(U_tot));
shading interp;
title(['Imaginary part of the total field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off

ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, abs(U_tot));
shading interp;
title(['Absolute value of the total field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off
