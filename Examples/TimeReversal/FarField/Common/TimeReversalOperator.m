% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute time reversal operator T = F*F = FF*
% in the far fielc and acoustic case
% Input
% Nb_angles        : nombre d'angles de discrétisation
% halpha           : pas de discrétisation des angles
% EmissionAngles  : vecteur des angles d'émission (rad)
% ReceptionAngles : vecteur des angles de réception (rad) (= EmissionAngles + pi)
% O [2 x N_scat]   : matrice des centres des obstacles
% a [1 x N_scat]   : vecteur des rayons des obstacles
% k                : fréquence
% Paramètres de sortie
% T : matrice de RT
% F : matrice de champ lointain
% M_modes : vecteurs des modes de Fourier
%

function [T, F, M_modes] = TimeReversalOperator(O, a, k, halpha, EmissionAngles, ReceptionAngles, varargin)

N_scat = length(a);
tolM = 10.^(-8);
%% Init Solution
Nb_anglesE = length(EmissionAngles);
Nb_anglesR = length(ReceptionAngles);
F = zeros(Nb_anglesR,Nb_anglesE);


%% Matrix of the system
nvar = length(varargin);
if (nvar == 0)
    error('Which type of condition ? varargin must be set (1=Dirichlet, 2=Neumann, 3=Penetrable)');
else
   if(isscalar(varargin{1}))
       if(varargin{1} == 1)
          Penetrable = 0;
          TypeOfOp = 6; %Dirichlet
          TypeOfForm = [1,0]; %Single layer potential
       elseif(varargin{1} == 2)
          Penetrable = 0;
          TypeOfOp =7; %Neumann
          TypeOfForm = [0,1]; % Double layer potential
       elseif(varargin{1} == 3)
          Penetrable = 1;
          omega = varargin{2};
          epsilon_0 = varargin{3};
          mu_0 = varargin{4};
          epsilon_m = varargin{5};
          mu_m = varargin{6};
          k_int = omega .* sqrt(epsilon_m .* mu_m);
          if(~isscalar(k_int) && length(k_int) < N_scat)
            error('vector k_int not large enough');              
          end
       else
            error('Unknown boundary condition');
       end
   else
        error('varargin{1} must be scalar (1=Dirichlet, 2=Neumann, 3= Penetrable)');
   end
end

%% Penetrable case
if(~Penetrable)
    % Fourier bases truncations
    M_modes = max(ones(1,N_scat),max([2^2*ones(1,N_scat); floor(k.*a + (1/(2*sqrt(2))*log(2*sqrt(2)*pi*k.*a/tolM)).^(2/3).*(k.*a).^(1/3) +1)]));
    %Matrix
    IntOp = PrecondDirichlet(O, a, M_modes, k);
    % Loop on incident waves
    for indice_angle = 1:Nb_anglesE
        %Angle of incidence
        beta_inc = EmissionAngles(indice_angle);
        %Right hand side
        [B] = PlaneWave(O, a, M_modes, k, beta_inc);
        %Density solution
        rho = IntOp\B;
        %Far Field
        F(:,indice_angle) = halpha*FarField(O, a, M_modes, k, ReceptionAngles, rho, TypeOfForm);
    end
%% Impenetrable case    
else
    %Building the integral operator
    M_modes_plus = max(10*ones(1,N_scat),max([2^2*ones(1,N_scat); floor(k.*a + (1/(2*sqrt(2))*log(2*sqrt(2)*pi*k.*a/tolM)).^(2/3).*(k.*a).^(1/3) +1)]));
    M_modes_minus = max(10*ones(1,N_scat),max([2^2*ones(1,N_scat); floor(k_int.*a + (1/(2*sqrt(2))*log(2*sqrt(2)*pi*k_int.*a/tolM)).^(2/3).*(k_int.*a).^(1/3) +1)]));
    % Truncation indices are chosen as the max between intern and extern values
    M_modes = max(M_modes_plus, M_modes_minus);
    sum_modes = sum(2*M_modes+1);
    %
    lambda_matrix = zeros(sum_modes,sum_modes);
    Sp = 0;
    for p = 1:N_scat
        Np = 2*M_modes(p)+1;
        lambda_matrix(Sp+1:Sp+Np, Sp+1:Sp+Np) = (mu_0/mu_m(p))*eye(Np);
        Sp = Sp +Np;
    end
    % Boundary integral operators
    Splus = SingleLayer(O, a, M_modes, k); 
    Sminus = IntegralOperator(O, a, M_modes, k_int, 2*eye(N_scat,N_scat));
    Nplus = DnSingleLayer(O, a, M_modes, k);
    Nminus = IntegralOperator(O, a, M_modes, k_int, 4*eye(N_scat,N_scat));
    Identity = eye(size(Nplus));
    % Matrix of the system
    A = [Splus, -Sminus; -0.5*Identity + Nplus, -lambda_matrix*(0.5*Identity + Nminus)];
    diagA = diag(1./(diag(A)));
    A = diagA*A;
    %cleaning memory
    clear Splus Sminus Nplus Nminus Identity lambda_matrix;
    % Loop on incident waves
    for indice_angle = 1:Nb_anglesE
        %Angle of incidence
        beta_inc = EmissionAngles(indice_angle);
        %Right hand side
        B1 = PlaneWave(O, a, M_modes, k, beta_inc);
        B2 = DnPlaneWave(O, a, M_modes, k, beta_inc);
        [B] = diagA*[B1;B2];
        clear B1 B2;
        %Density solution
        rho = A\B;
        rho_plus = rho(1:sum_modes);
        %Far Field
        F(:,indice_angle) = halpha*FarField(O, a, M_modes, k, ReceptionAngles, rho_plus, [1,0]);
    end    
end
%% Time reversal operator
T = F'*F;

