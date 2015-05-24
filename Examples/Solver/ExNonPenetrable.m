% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Contributed by X. Antoine and B. Thierry
%
% Example of a solution of a Dirichlet or Neumann problem (non penetrable)
% using the DirichletSolver (resp. NeumannSolver) of mu-diff.
% 
disp('-----------------------');
disp('    New simulation     ');
disp('-----------------------');
clear all;
close all;
ind_fig = 1;

%% Chose the boundary condition
TypeOfProblem = 'Dirichlet';
%TypeOfProblem = 'Neumann';

disp(['Boundary condition: ', TypeOfProblem]);
%% GRID on which the fiels are computed
XXmin = -10; XXmax = 10;
YYmin = -10; YYmax = 10;
lc = min(2*pi/15, 1/10); % characteristic length

%% -----------------------------------------
%%            Pre Processing
%% -----------------------------------------
%% Incident Wave (chose one)
TypeOfWave = 'PlaneWave';
%TypeOfWave = 'PointSource';
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
%% Wavenumber
k =2*pi;
disp(['Wavenumber: ', num2str(k)]);
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
%% Drawing the circular obstacles
ind_fig = ind_fig +1;
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2);
axis equal;
xlabel('x'); ylabel('y');
title('Obstacles');

%% -----------------------------------------
%% Solving the problem
%% -----------------------------------------
if(strcmp(TypeOfProblem, 'Dirichlet'))
    Solution = DirichletSolver(O, a, k, TypeOfWave, ParamWave);
else
    Solution = NeumannSolver(O, a, k, TypeOfWave, ParamWave);
end    
%% -----------------------------------------
%%            Post Processing
%% -----------------------------------------
%% Radar Cross Section
theta_RCS = 0:360;
theta_RCS_rad = theta_RCS*2*pi/360;

if(strcmp(TypeOfProblem, 'Dirichlet'))
    R = DirichletRCS(Solution, theta_RCS_rad);
else
    R = NeumannRCS(Solution, theta_RCS_rad);
end
ind_fig = ind_fig +1;
figure(ind_fig);
plot(theta_RCS_rad, R, 'k-');
hold on
title('Radar cross section');
xlabel('Angle of reception (degree)');
ylabel('Radar Cross Section (dB)');
axis tight;

%% Field on grid (based on an integral equation representation)
%SEE BEGINING FOR XXmin, XXmax, YYmin and YYmax values
XX = [XXmin:lc:XXmax];
YY = [YYmin:lc:YYmax];
[X,Y] = meshgrid(XX,YY);

if(strcmp(TypeOfProblem, 'Dirichlet'))
    [U_tot, U] = DirichletNearField(Solution, X, Y);
else
    [U_tot, U] = NeumannNearField(Solution, X, Y);
end    
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
