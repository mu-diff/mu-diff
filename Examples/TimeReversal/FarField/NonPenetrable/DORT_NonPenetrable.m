% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
%-------------------------
% DORT acoustic in 2D.
% Contributed by X. Antoine, K. Ramdani and B. Thierry
% Related articles:
%1) C. Hazard, K. Ramdani.
%   "Selective acoustic focusing using time-harmonic reversal mirrors."
%   SIAM J. Appl. Math., (2004)
%
%2) B. Thierry
%   "Analyse et Simulations Numériques du Retournement Temporel et de la 
%    Diffraction Multiple"
%   PhD Thesis, University of Nancy, France, 2011 (in french).
%

%%
disp('------- New simulation ------');
clear all;
close all;

ind_fig = 1;
%Wavelength and wavenumber
lambda =1;
k = 2*pi/lambda;
CONDITION = 1; % Dirichlet =1, Neumann = 2

%% Geometry
O = [0, 10, -10;20, -10, -20];
a =[0.02, 0.01, 0.005];

N_scat = size(O,2);

%sorting in descending order
[a, IndexList] = sort(a,2, 'descend');
O = O(:,IndexList);

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
T = TimeReversalOperator(O, a, k, halpha, EmissionAngles, ReceptionAngles, CONDITION);
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
%    print('-dpng','Eigenvalue');
end

%% -----------------------------
% Herglotz waves of eigenvectors
% ------------------------------
%Mesh
step_d_x = 1/10;
step_d_y = 1/10;
xmin = -val_max; xmax = val_max;
ymin = -val_max; ymax = val_max;
[X,Y] = meshgrid(xmin:step_d_x:xmax,ymin:step_d_y:ymax);
[nby,nbx] = size(X);
%% Eigvenvectors
%Number of Herglotz wave to compute (time-consuming !)
Nb_vect = 3;
disp(['Launching the computation of the ',num2str(Nb_vect),' Herglotz waves']);

cell_vectT = cell(Nb_vect,1);

for cpt_vect = 1:(Nb_vect)
    disp(['Herglotz wave numbered ',num2str(cpt_vect),' of ',num2str(Nb_vect)]);    
    V = vect(:,end+ 1- cpt_vect);
    cell_vectT{cpt_vect} = HerglotzWave(V, k, halpha, EmissionAngles, X,Y);
end

disp([num2str(Nb_vect),' Herglotz waves: done']);

%% Display
for cpt_vect = 1:Nb_vect
    Z =  cell_vectT{cpt_vect};
    ind_fig = ind_fig+1; figure(ind_fig);
    pcolor(X,Y,abs(Z)); %view(2); 
    colorbar;
    shading interp;
    title(['Herglotz wave associated to eigenvector n°',NUM2STR(cpt_vect)]) 
    xlabel('x_1'); ylabel('x_2');
    axis([xmin, xmax, ymin,ymax]);
    hold on 
%    plot(O(1, 1:N_scat), O(2, 1:N_scat), 'wo');
    hold off
%    print('-dpng',['Herglotz',num2str(cpt_vect)]);
end

%% Some values

for cpt = 1: max(Nb_vect, N_scat+1)
    disp(['Eigenvalue ',num2str(cpt),' : ', num2str(eigenval(end - cpt+1))]);
end

%% Max of Herglotz waves
O_maxZ = zeros(2,Nb_vect);
for cpt = 1:Nb_vect
    Z = abs(cell_vectT{cpt});
    maxZ = max(max(Z));
    [maxZ_x,maxZ_y] = find(Z == maxZ);
    O_maxZ(:,cpt) = [X(maxZ_x,maxZ_y);Y(maxZ_x,maxZ_y)];
end

disp('Localisation of the maximal values of the Herglotz waves...');
for cpt = 1:Nb_vect
    disp(['Herglotz wave numbered ',num2str(cpt),': (', num2str(O_maxZ(1,cpt)),', ',num2str(O_maxZ(2,cpt)),')']);
end

%%
disp('End of simulation');

