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
%load('data.mat');
%load('dataN1.mat');
%load('data10.mat');
%load('dataGMRES.mat');
load('dataGMRESM4.mat')
%load('dataDirect.mat');


ind_fig = 1;

%% Drawing the circular obstacles
PlotCircles(O, a, ind_fig, 'Color', 'k', 'LineWidth', 2);
%axis([x_min, x_max, y_min, y_max]);
ind_fig = ind_fig +1;
xlabel('x'); ylabel('y');
title('Obstacles');
axis equal;
%% Radar Cross Section (RCS)
ind_fig = ind_fig +1; figure(ind_fig);
plot(theta_RCS, R, ColorBIE);
title('Radar Cross Section');
xlabel('Angle of reception (degree)');
ylabel('Radar Cross Section (dB)');
axis tight;

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
