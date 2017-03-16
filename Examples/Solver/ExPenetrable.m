% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
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
% Example of a solution of a penetrable problem 
% using the PenetrableSolver of mu-diff.
% 

disp('-----------------------');
disp('    New simulation     ');
disp('-----------------------');

clear all;
close all;
ind_fig = 1;

%% -----------------------------------------
%%            Pre Processing
%% -----------------------------------------
%% Incident Wave (chose one)
%TypeOfWave = 'PlaneWave';
TypeOfWave = 'PointSource';
switch TypeOfWave
    case 'PlaneWave',
        % Incidence angle
        ParamWave = pi;
    case 'PointSource',
        %coordinated of the point source
        ParamWave = [-10; 0];
    otherwise,
        error('Unknown type of wave');
end
disp(['Incident wave: ', TypeOfWave]);

%% GRID on which the fiels are computed
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
contrast = [2.2, 2.2];
k = 1;  %exterior wavenumber
k_int = k.*contrast;

disp(['Wavenumber (exterior): ', num2str(k)]);
str_disp = num2str(k_int(1));
for i=2:N_scat
    str_disp = [str_disp, ', ', num2str(k_int(i))];
end
disp(['Wavenumber (interior): ', str_disp]);

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
%% Solving the problem
%% -----------------------------------------
Solution = PenetrableSolver(O, a, k, k_int, TypeOfWave, ParamWave);

%% Radar Cross Section (RCS)
theta_RCS = 0:360;
theta_RCS_rad = theta_RCS*2*pi/360;

R = PenetrableRCS(Solution, theta_RCS_rad);

ind_fig = ind_fig +1; figure(ind_fig);
plot(theta_RCS, R, 'k-');
title('Radar Cross Section');
xlabel('Angle of reception (degree)');
ylabel('Radar Cross Section (dB)');
axis tight;

%% Scattered field on grid
%SEE BEGINING FOR XXmin, XXmax, YYmin and YYmax values
XX = [XXmin:lc:XXmax];
YY = [YYmin:lc:YYmax];
[X,Y] = meshgrid(XX,YY);

[U_tot, U] = PenetrableNearField(Solution, X, Y);


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
