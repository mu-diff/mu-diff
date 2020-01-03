% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
%-------------------------
% DORT dielectric in 2D.
% Contributed by C. Burkard, K. Ramdani and B. Thierry
% Related article:
%   C. Burkard, A. Minut, K. Ramdani, 
%   "Far field model for time reversal and application to selective focusing on 
%    small dielectric inhomogeneities."
%   Inverse Problems and Imaging, 2013.
%-------------------------

clear all;
close all;
format long;

disp('------- Nouveau Test ------');

%% Setup 
ind_fig = 1;

%domain......
xmin = -5;
xmax = 5;
ymin = 0;
ymax = 20;

%Wavenumber and physical parameters
epsilon_0    = 1;
mu_0         = 1;
omega       = 2*pi;
k           = 2*pi*sqrt(epsilon_0 * mu_0); 
lambda      = 2*pi/k;

%interior constants and wave numbers ..............
epsilon_m   = [10*epsilon_0 , 5*epsilon_0 ];     %, lambda + 0.5, lambda + 3, lambda + 2];
mu_m        = [mu_0    , mu_0 ];      %, lambda + 2, lambda + 2, lambda + 1];
k_m         = omega .* sqrt(epsilon_m .* mu_m);



%% Geometry (scatterers)
O = [0 , 2 ;10 , 2];%,0,4, 8; 6, 17,11, 2];
a = [0.03  ,  0.03];
%sorting in descending order
[a, IndexList] = sort(a,2, 'descend');
O = O(:,IndexList);
N_scat = size(O,2);


%% Plot of the Scatterers
ind_fig = ind_fig +1; figure(ind_fig);
% disks are too small to be displayed with PlotCircles, only their center is
% thus drawn
LineWidth = [N_scat:-1:1];
hold on
for p=1:N_scat
    plot(O(1,p),O(2,p),'ko','LineWidth', LineWidth(p));    
end
hold off
xlabel('x_1'); ylabel('x_2');
title('Obstacles (center of)');
%Usefull later for Herglotz waves
AxisValues = axis;
val_max = max(abs(AxisValues))+2;
axis([-val_max,val_max,-val_max,val_max])

%% Test for order of error in DORT O(1/kd).........
distScatterers = zeros(N_scat,N_scat);
for q= 1:N_scat
    x0 = O(1,q);
    y0 = O(2,q);
    for p = 1:N_scat
        x = O(1,p);
        y = O(2,p);
        if q==p
            distScatterers(q,p) = 100;
        else
            distScatterers(q,p) = norm([x-x0  y-y0]);
        end
    end
end
distmin = min(min(distScatterers));
testpar = 1/(k * distmin);
disp(['Error order  O(1/(kd))  =    ' , num2str(testpar)]);

%% Building the TRM

%Reception angles (start - end)
ReceptionAnglesStart = 0;
ReceptionAnglesEnd = 2*pi;

% TRM size (can be reduce to get an open TRM)
TRM_size = abs(ReceptionAnglesEnd  - ReceptionAnglesStart);
TRM_size_deg = TRM_size *360/2/pi; %degree version

%Number of discretization points of the TRM
Nb_angles = 720;
%discretization step of [0:2pi]
halpha = TRM_size/(Nb_angles);
halpha_deg = TRM_size_deg/(Nb_angles);

%Reception angles vector of the TRM
if (mod(ReceptionAnglesStart- ReceptionAnglesStart, 2*pi) == 0)
    %closed TRM then delete doubled angle
    ReceptionAngles = ReceptionAnglesStart:halpha:(ReceptionAnglesEnd-halpha);
else    
    %open TRM
    ReceptionAngles = ReceptionAnglesStart:halpha:ReceptionAnglesEnd;
    Nb_angles = Nb_angles +1;
end
ReceptionAngles_deg = ReceptionAngles * (360/2/pi);

%Emission angles of the TRM
EmissionAngles = ReceptionAngles + pi;
EmissionAngles_deg = ReceptionAngles_deg +180;

%Possible dephasing
phi = 0; phi_deg = phi*360/2/pi;

ReceptionAngles_deg = ReceptionAngles_deg + phi_deg;
EmissionAngles_deg = EmissionAngles_deg + phi_deg;
ReceptionAngles = ReceptionAngles + phi;
EmissionAngles = EmissionAngles + phi;

%% ----------------------
% Time reversal operator
% -----------------------                       
T = TimeReversalOperator(O, a, k, halpha, EmissionAngles, ReceptionAngles, 3, omega, epsilon_0, mu_0, epsilon_m, mu_m);
% Eigenvalues
[vect,val] = eig(T);
disp('Time reversal matrix and eigen values: done');

%% ------------------------
% Study of the eigenvalues
% -------------------------
eigenval = sort(abs(diag(val)));
%The N_scat^{th} will be plotted in red
val_N_scat = eigenval(Nb_angles + 1 - N_scat);

% Normal scale
ind_fig = ind_fig +1; figure(ind_fig);
hold on
plot([Nb_angles:-1:1], eigenval,'*')
plot(N_scat,val_N_scat,'r*');
hold off
xlabel('Index of the eigenvalues');
ylabel('Absolute value of the eigenvalues');
title('Absolute value of the eigenvalues of the time reversal operator');

% Log scale
ind_fig = ind_fig +1; figure(ind_fig); 
hold on
plot([Nb_angles:-1:1], log10(eigenval),'*')
plot(N_scat,log10(val_N_scat),'r*');
hold off 
xlabel('Index of the eigenvalues');
ylabel('Absolute value of the eigenvalues (LOG_{10})');
title('Absolute value of the eigenvalues of the time reversal operator');

%Zoom on the first 40 eigenvalues
if Nb_angles > 40
    val_vp_zoom = 40;
    tab_zoom_vp = 1:val_vp_zoom ;

    ind_fig = ind_fig +1; figure(ind_fig);
    hold on
    plot([val_vp_zoom:-1:1], log10(eigenval(end+1-val_vp_zoom:end)),'*')
    plot(N_scat,log10(val_N_scat),'r*');
    hold off
    xlabel('Index of the eigenvalues');
    ylabel('Absolute value of the eigenvalues (LOG_{10})');
    title('Absolute value of the eigenvalues of the time reversal operator');
end

%% Significant eigenvalues
v_sig = 0;v_Sig = 0;
N_vp = length(eigenval);
for ell=1:N_vp
    if (eigenval(end+1-ell) > 10^(-6) )
        v_Sig(ell) = ell;
    end
end

Nb_vect = N_scat*(2 + 1);

ind_fig= ind_fig+1; figure(ind_fig);
plot(1:N_vp, log10(eigenval(end:-1:1)), 'b*');
hold on
plot(1:Nb_vect, log10(eigenval(end:-1:end+1-Nb_vect)), 'r.');
hold off;

disp(['Number of significant Eigenvalues: ' num2str(length(v_Sig))]);


%% Mesh for Herglotz waves
step_x = lambda/30;
step_y = lambda/30;
[X,Y]   = meshgrid(xmin:step_x:xmax,ymin:step_y:ymax);
[nby,nbx] = size(X);
%% Herglotz waves
disp(['Launching the computation of the ',num2str(Nb_vect),' Herglotz waves']);
cell_vectT = cell(Nb_vect,1);
for cpt_vect = 1:(Nb_vect)
    disp(['Herglotz wave numbered ',num2str(cpt_vect),' of ',num2str(Nb_vect)]);    
    V = vect(:,end+ 1- cpt_vect);
    cell_vectT{cpt_vect} = HerglotzWave(V, k, halpha, EmissionAngles, X,Y);
end

disp([num2str(Nb_vect),' Herglotz waves: done']);

%% Draw Herglotz Waves
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2]);
drawnow;

for cpt_vect = 1:(Nb_vect)
    Z = abs(cell_vectT{cpt_vect});
    subplot(1,Nb_vect, cpt_vect);
    pcolor(X,Y,Z);
    hold on 
    plot(O(1, 1:N_scat), O(2, 1:N_scat), 'wo');
    colorbar;
    title(['Herglorz wave number ',num2str(cpt_vect)]); 
    xlabel('x_1'); ylabel('x_2');
    axis tight;
    shading interp;
end




