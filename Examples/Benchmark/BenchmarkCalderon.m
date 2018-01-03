% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
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
%% BIE: first (KIND=1) or second (KIND=2) Kind ?
KIND = 1;

if(KIND == 1)
    disp('Integral equation of the first kind')
elseif(KIND==2)
    disp('Integral equation of the second kind')
else
    error('KIND must be equal to 1 or 2');
end
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
% Helmholtz representation:
% u = -/+ \int_{\Gamma} G^{+/-}(x,y) u(y) d\Gamma(y)
%          -/+ \int_{\Gamma} dn G^{+/-}(x,y) dn u(y) d\Gamma(y)
% Calderon operator
%% Boundary integral operators
Cext = CalderonProjector(O, a, M_modes, k);
Cint = CalderonProjector(O, a, M_modes, k_int, eye(N_scat));
Cext = Cext;
Cint = Cint - eye(size(Cint));
switch KIND
    case 1,
        A = Cext - Cint;
    case 2,
        A = Cext + Cint;
end
clear Cext Cint;
%% Incident plane wave (right hand side)
switch IncidentWave
    case 'PlaneWave',
        Uinc = PlaneWave(O, a, M_modes, k, beta_inc);
        DnUinc = DnPlaneWave(O, a, M_modes, k, beta_inc);
    case 'PointSource',
        Uinc = PointSource(O, a, M_modes, k, XS);
        DnUinc = DnPointSource(O, a, M_modes, k, XS);
end
%Reshape
B = zeros(size([Uinc;DnUinc]));
Sp_uinc = 0;
Sp_B = 0;
for p =1:N_scat
    Np = M_modes(p);
    MNp = [-Np:Np];
    B(Sp_B + MNp + 1 + Np,1) = Uinc(Sp_uinc + MNp + 1 + Np,1);
    Sp_B = Sp_B + 2*Np+1;
    B(Sp_B + MNp + 1 + Np,1) = DnUinc(Sp_uinc + MNp + 1 + Np,1);
    Sp_B = Sp_B + 2*Np+1;
    Sp_uinc = Sp_uinc + 2*Np+1;
end
B = -B;
clear Uinc DnUinc;
%% Solving linear system
AllTraces = A\B;
%% Extracting the traces and normal derivative traces
gammaD = zeros(sum_modes,1);
gammaN = zeros(sum_modes,1);
Sp_gamma = 0;
Sp_Matrix = 0;
for p =1:N_scat
    Np = M_modes(p);
    MNp = [-Np:Np];
    gammaD(Sp_gamma + MNp + 1 + Np,1) = AllTraces(Sp_Matrix + MNp + 1 + Np);
    Sp_Matrix = Sp_Matrix + 2*Np+1;
    gammaN(Sp_gamma + MNp + 1 + Np,1) = AllTraces(Sp_Matrix + MNp + 1 + Np);
    Sp_Matrix = Sp_Matrix + 2*Np+1;
    Sp_gamma = Sp_gamma + 2*Np+1;
end
%% Radar Cross Section (RCS)
theta_RCS = 0:360;
theta_RCS_rad = theta_RCS*2*pi/360;
R = RCS(O, a, M_modes, k, theta_RCS_rad, [gammaN, gammaD], [-1,-1]);
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

Ue = ExternalPotential(X, Y, O, a, M_modes, k, [gammaN,gammaD], [-1,-1]);
Ui = InternalPotential(X, Y, O, a, M_modes, k_int, [gammaN,gammaD], [1,1], 'OnBoundary', 1);

%% Incident wave on grid (outside obstacles)
MeshOutside = MaskMatrixObstacles(X, Y, O, a) == 0;
UincOnMesh = zeros(size(X));
switch IncidentWave
    case 'PlaneWave',
        XYBeta = X*cos(beta_inc)+ Y*sin(beta_inc);
        UincOnMesh(MeshOutside) = exp(1i*k*XYBeta);
    case 'PointSource',
        ABS_XY = sqrt((X-XS(1)).^2+ (Y-XS(2)).^2);
        UincOnMesh(MeshOutside) = 1i/4*besselh(0, 1, k*ABS_XY(MeshOutside));
end
%% Total field
U_tot = Ue + Ui + UincOnMesh;
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
