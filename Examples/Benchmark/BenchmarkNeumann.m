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
% incident wave by a collection of sound-hard (Neumann) circular obstacles.
% The Benchmark file solves the problem with four different integral equations
% (EFIE, MFIE, CFIE and Brakhage-Werner), with both dense and sparse
% versions. The single-scattering preconditioned integral equation is
% moreover also solved, for the dense and sparse versions.
% --------------------------------------------------------------------
% The Radar Cross Section (RCS) is also computed for every chosen integral
% equations, and the scattered field is computed on a mesh for only one
% chosen integral formulation (can be time-consuming for large domains).
% 
disp('-----------------------');
disp('    New simulation     ');
disp('-----------------------');
clear all;
close all;
%% Some parameters
ind_fig = 1;
%which integral equations are solved ?
%(only use as a convenient way)
int_equation = {};

int_equation=[int_equation, 'EFIE'];
int_equation=[int_equation, 'MFIE'];
int_equation=[int_equation, 'CFIE'];
int_equation=[int_equation, 'BW'];
int_equation=[int_equation, 'SpEFIE'];
int_equation=[int_equation, 'SpMFIE'];
int_equation=[int_equation, 'SpCFIE'];
int_equation=[int_equation, 'SpBW'];
int_equation=[int_equation, 'Precond'];
int_equation=[int_equation, 'SpPrecond'];
%% colors
ColorEFIE = 'k-o';
ColorMFIE = 'r-o';
ColorCFIE = 'm-o';
ColorBW = 'b-o';
ColorPrecond = 'c-o';
ColorSpEFIE = 'k-.';
ColorSpMFIE= 'r-.';
ColorSpCFIE= 'm-.';
ColorSpBW = 'b-.';
ColorSpPrecond = 'c-.';

%% GRID on which the fiels are computed
XXmin = -10; XXmax = 10;
YYmin = -10; YYmax = 10;
%lc = min(lambda/15, 1/10); % characteristic length
lc = min(2*pi/15, 1/10); % characteristic length

%% -----------------------------------------
%%            Pre Processing
%% -----------------------------------------
%% Incident Wave (chose one)
IncidentWave = 'PlaneWave';
%IncidentWave = 'PointSource';
% Incidence angle
beta_inc = pi;
%For point source
XS = [-10; 0];
%% Wavenumber
k =2*pi;
%Wavelength
lambda = 2*pi/k;

%% Geometry
GEOMETRY = 'manual';
%GEOMETRY = 'random';
%GEOMETRY = 'rectangular';
%GEOMETRY = 'triangular';

switch GEOMETRY
    case 'manual',
        O = [-2.5, 2.5    ; -2, 1];
        a = [1.5, 0.75];
    case 'random',
        x_min = -5; x_max = 5;
        y_min = -5; y_max = 5;
        distance_min = 0.01;
        a_min = 0.9;
        a_max = 1.1;
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
axis equal;
ind_fig = ind_fig +1;
xlabel('x'); ylabel('y');
title('Obstacles');

%% Fourier series: truncation indices
M_modes = FourierTruncation(a, k, 'Min', 1);

%% GMRES parameters
RESTART = [];
TOL = 10^(-10);
MAXIT = sum(2*M_modes+1);
%% Incident plane wave (right hand side)
switch IncidentWave
    case 'PlaneWave'
        Uinc = PlaneWave(O, a, M_modes, k, beta_inc);
        DnUinc = DnPlaneWave(O, a, M_modes, k, beta_inc);
    case 'PointSource'
        Uinc = PointSource(O, a, M_modes, k, XS);
        DnUinc = DnPointSource(O, a, M_modes, k, XS);
end

%% -----------------------------------------
%%            Integral formulations
%% -----------------------------------------
%% EFIE
A_EFIE = DnDoubleLayer(O, a, M_modes, k);
[rho_EFIE,FLAG_EFIE,RELRES_EFIE,ITER_EFIE,RESVEC_EFIE] = gmres(A_EFIE, DnUinc, RESTART, TOL, MAXIT, [], []);
%% SPARSE EFIE
SpA_EFIE = SpDnDoubleLayer(O, a, M_modes, k);
[rho_SpEFIE,FLAG_SpEFIE,RELRES_SpEFIE,ITER_SpEFIE,RESVEC_SpEFIE] = gmres(@(X)SpMatVec(X, M_modes, SpA_EFIE), DnUinc, RESTART, TOL, MAXIT, [], []);

%% MFIE
A_MFIE = IntegralOperator(O, a, M_modes, k, {'I','M'}, [0.5, 1]);
[rho_MFIE,FLAG_MFIE,RELRES_MFIE,ITER_MFIE,RESVEC_MFIE] = gmres(A_MFIE, Uinc, RESTART, TOL, MAXIT, [], []);
%% SPARSE MFIE
SpM = SpDoubleLayer(O, a, M_modes, k);
SpI = SpIdentity(O, a, M_modes);
%[rho_SpMFIE,FLAG_SpMFIE,RELRES_SpMFIE,ITER_SpMFIE,RESVEC_SpMFIE] = gmres(@(X)SpMatVec(X,M_modes, {SpI, SpA_MFIE}, [0.5,1]), DnUinc, RESTART, TOL, MAXIT, [], []);
%Or compute directly the sparse matrix:
SpA_MFIE = SpAddIdentity(SpM, 0.5, M_modes);
[rho_SpMFIE,FLAG_SpMFIE,RELRES_SpMFIE,ITER_SpMFIE,RESVEC_SpMFIE] = gmres(@(X)SpMatVec(X, M_modes, SpA_MFIE), Uinc, RESTART, TOL, MAXIT, [], []);

%% CFIE
alpha = 0.5;
eta = -1i*k;
A_CFIE = (1-alpha) * A_EFIE + eta*alpha* A_MFIE;
B_CFIE = alpha * eta * Uinc + (1-alpha)* DnUinc;
[rho_CFIE,FLAG_CFIE,RELRES_CFIE,ITER_CFIE,RESVEC_CFIE] = gmres(A_CFIE, B_CFIE, RESTART, TOL, MAXIT, [], []);
%% SPARSE CFIE
[rho_SpCFIE, FLAG_SpCFIE, RELRES_SpCFIE, ITER_SpCFIE, RESVEC_SpCFIE] = gmres(@(X)SpMatVec(X, M_modes,{SpA_EFIE, SpA_MFIE}, [1-alpha, alpha*eta]), B_CFIE, RESTART, TOL, MAXIT, [], []);

%% Brackage Werner
eta_BW = 1i/k;
A_BW = IntegralOperator(O, a, M_modes, k, {'I', 'N', 'D'}, [0.5, -1, -eta_BW]);
[psi_BW,FLAG_BW,RELRES_BW,ITER_BW,RESVEC_BW] = gmres (A_BW, DnUinc, RESTART,TOL,MAXIT);
%% SPARSE Brackage Werner
SpN = SpDnSingleLayer(O, a, M_modes, k);
[psi_SpBW,FLAG_SpBW,RELRES_SpBW,ITER_SpBW,RESVEC_SpBW] = gmres(@(X)SpMatVec(X, M_modes,{SpI, SpN, SpA_EFIE}, [0.5, -1, -eta_BW]), DnUinc, RESTART, TOL, MAXIT, [], []);

%% Single scattering preconditioned integral equation (single-layer formulation)
A_Precond = IntegralOperator(O, a, M_modes, k, {'Dprec'});
switch IncidentWave
    case 'PlaneWave'
        DnUincPrecond = DnPlaneWavePrecond(O, a, M_modes, k, beta_inc);
    case 'PointSource'
        DnUincPrecond = diag(1./diag(A_EFIE))*DnUinc;
end

[rho_Precond,FLAG_Precond,RELRES_Precond,ITER_Precond,RESVEC_Precond] = gmres(A_Precond, DnUincPrecond, RESTART, TOL, MAXIT, [], []);
%% SPARSE Precond
SpPrecond = SpIntegralOperator(O, a, M_modes, k, {'Dprec'});
[rho_SpPrecond,FLAG_SpPrecond,RELRES_SpPrecond,ITER_SpPrecond,RESVEC_SpPrecond] = gmres(@(X)SpMatVec(X, M_modes,SpPrecond), DnUincPrecond, RESTART, TOL, MAXIT, [], []);

%% -----------------------------------------
%%            Post Processing
%% -----------------------------------------
%% History of convergence of the GMRES
figure(ind_fig); ind_fig = ind_fig +1;
Legende = {};
hold on
for cpt_ie = 1:length(int_equation)
    switch int_equation{cpt_ie}
        case 'EFIE',
            resv = RESVEC_EFIE; iterv = ITER_EFIE;
            res = log10(RESVEC_EFIE./max(RESVEC_EFIE));
            col = ColorEFIE; leg = 'EFIE';
        case 'MFIE',
            resv = RESVEC_MFIE; iterv = ITER_MFIE;
            res = log10(RESVEC_MFIE./max(RESVEC_MFIE));
            col = ColorMFIE; leg = 'MFIE';
        case 'CFIE',
            resv = RESVEC_CFIE; iterv = ITER_CFIE;
            res = log10(RESVEC_CFIE./max(RESVEC_CFIE));
            col = ColorCFIE; leg = 'CFIE';
        case 'BW',
            resv = RESVEC_BW; iterv = ITER_BW;
            res = log10(RESVEC_BW./max(RESVEC_BW));
            col = ColorBW; leg = 'BW';
        case 'SpEFIE',
            resv = RESVEC_SpEFIE; iterv = ITER_SpEFIE;
            res = log10(RESVEC_SpEFIE./max(RESVEC_SpEFIE));
            col = ColorSpEFIE; leg = 'Sparse EFIE';
        case 'SpMFIE',
            resv = RESVEC_SpMFIE; iterv = ITER_SpMFIE;
            res = log10(RESVEC_SpMFIE./max(RESVEC_SpMFIE));
            col = ColorSpMFIE; leg = 'Sparse MFIE';
        case 'SpCFIE',
            resv = RESVEC_SpCFIE; iterv = ITER_SpCFIE;
            res = log10(RESVEC_SpCFIE./max(RESVEC_SpCFIE));
            col = ColorSpCFIE; leg = 'Sparse CFIE';
        case 'SpBW',
            resv = RESVEC_SpBW; iterv = ITER_SpBW;
            res = log10(RESVEC_SpBW./max(RESVEC_SpBW));
            col = ColorSpBW; leg = 'Sparse BW';
        case 'Precond',
            resv = RESVEC_Precond; iterv = ITER_Precond;
            res = log10(RESVEC_Precond./max(RESVEC_Precond));
            col = ColorPrecond; leg = 'Precond';
        case 'SpPrecond',
            resv = RESVEC_SpPrecond; iterv = ITER_SpPrecond;
            res = log10(RESVEC_SpPrecond./max(RESVEC_SpPrecond));
            col = ColorSpPrecond; leg = 'Sparse Precond';
    end
    if(RESTART >= 1 )
        plot(1:length(resv), res, col); Legende = [Legende; leg];
    else
        plot(1:length(resv), res, col); Legende = [Legende; leg];
    end
end

if( RESTART >= 1 )
    title(['History of convergence of the GMRES, k = ',num2str(k), ', restart = ', num2str(RESTART)])
else
    title(['History of convergence of the GMRES, k = ',num2str(k), ', without restart'])
end
legend(Legende);
xlabel('Iteration');
ylabel('Residu (log)');

%% Radar cross section
theta_RCS = 0:360;
theta_RCS_rad = theta_RCS*2*pi/360;

figure(ind_fig);
ind_fig = ind_fig +1;
Legende = {};
hold on
for cpt_ie = 1:length(int_equation)
    switch int_equation{cpt_ie}
        case 'EFIE',
            RCS_EFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_EFIE, [0,1]);
            plot(theta_RCS, RCS_EFIE, ColorEFIE); Legende = [Legende;'EFIE'];
        case 'MFIE',
            RCS_MFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_MFIE, [0,1]);
            plot(theta_RCS, RCS_MFIE, ColorMFIE); Legende = [Legende;'MFIE'];
        case 'CFIE',
            RCS_CFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_CFIE, [0,1]);
            plot(theta_RCS, RCS_CFIE, ColorCFIE); Legende = [Legende;'CFIE'];
        case 'BW',
            RCS_BW = RCS(O, a, M_modes, k, theta_RCS_rad, psi_BW, [-1,-eta_BW]);
            plot(theta_RCS, RCS_BW, ColorBW); Legende = [Legende;'BW'];
        case 'SpEFIE',
            RCS_SpEFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_SpEFIE, [0,1]);
            plot(theta_RCS, RCS_SpEFIE, ColorSpEFIE); Legende = [Legende;'SpEFIE'];
        case 'SpMFIE',
            RCS_SpMFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_SpMFIE, [0,1]);
            plot(theta_RCS, RCS_SpMFIE, ColorSpMFIE); Legende = [Legende;'SpMFIE'];
        case 'SpCFIE',
            RCS_SpCFIE = RCS(O, a, M_modes, k, theta_RCS_rad, rho_SpCFIE, [0,1]);
            plot(theta_RCS, RCS_SpCFIE, ColorSpCFIE); Legende = [Legende;'SpCFIE'];
        case 'SpBW',
            RCS_SpBW = RCS(O, a, M_modes, k, theta_RCS_rad, psi_SpBW, [-1,-eta_BW]);
            plot(theta_RCS, RCS_SpBW, ColorSpBW); Legende = [Legende;'SpBW'];
        case 'Precond',
            RCS_Precond = RCS(O, a, M_modes, k, theta_RCS_rad, rho_Precond, [0,1]);
            plot(theta_RCS, RCS_Precond, ColorPrecond); Legende = [Legende;'Precond'];
        case 'SpPrecond',
            RCS_SpPrecond = RCS(O, a, M_modes, k, theta_RCS_rad, rho_SpPrecond, [0,1]);
            plot(theta_RCS, RCS_SpPrecond, ColorSpPrecond); Legende = [Legende;'Sparse Precond'];
    end
end
legend(Legende)
title('Radar cross section');
xlabel('Angle of reception (degree)');
ylabel('Radar Cross Section (dB)');
axis tight;

%% Map of the waves (based on one integral equation)
%%Grid
XX = [XXmin:lc:XXmax];
YY = [YYmin:lc:YYmax];
[X,Y] = meshgrid(XX,YY);
Matrice_Obstacles = MaskMatrixObstacles(X, Y, O, a);
Matrice_Not_Obstacles = (Matrice_Obstacles == 0);

%%Scattered field
U = ExternalPotential(X, Y, O, a, M_modes, k, rho_SpPrecond, [0,1]);
%% Incident field
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

ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, abs(U));
shading interp;
title(['Absolute value of the scattered field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off

ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, real(U));
shading interp;
title(['Real part of the scattered field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off

ind_fig = ind_fig +1; figure(ind_fig);
hold on
surf(X,Y, imag(U));
shading interp;
title(['Imaginary part of the scattered field']);
xlabel('x'); ylabel('y');
view(2); colorbar;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2, 'zdata', max(max(abs(U_tot))));
set(gcf,'Renderer','Zbuffer');
hold off
