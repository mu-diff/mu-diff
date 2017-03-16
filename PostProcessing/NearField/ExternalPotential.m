% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute a linear combination of the single- and double-layer potential for 
% the SAME potential, on a meshgrid (X, Y) and for the multiple scattering
% problem by circular obstacles in the Fourier bases.
% -------------------------------------------------------------------------
% The computation is done OUTSIDE the obstacles and the returned variable U
% inside the obstacles is set to 0.
% -------------------------------------------------------------------------
% Single-Layer potential L_p rho_p associated to obstacle p, for a density rho_p
% in H^{-1/2}(Gamma_p) is given by
%    L_p rho_p(x) = \int_{Gamma_p} G(x,y) rho_p(y) dy, for all x in R^2\ Omega_p
% Double-Layer potential M_p lambda_p for a density lambda_p in H^{1/2}(Gamma_p):
%    M lambda(x) = - \int_{Gamma_p} dn_y G(x,y) lambda_p(y) dy, for all x in R^2 \ Omega_p
%
% Omega_p: obstacle numbered p.
% Gamma_p: boundary of Omega_p.
% -------------------------------------------------------------------------
% The resulting potential is a linear combination of the single- and
% double-layer potentials, but associated to the SAME potential Density_p:
%    U = sum_{p=1}^{N_scat} (alpha_p L_p + beta_p M_p) Density_p
% -------------------------------------------------------------------------
% U = ExternalPotential(X, Y, N_scat, O, a, k, M_modes, Density, TypeOfOperator)
%
% Input Arguments (N_scat = number of obstacles):
% ----------------
% X              [N_Y, N_X]   : X matrix of the meshgrid (N_Y = Nb. points
%                               in Y-direction and N_X = Nb. points in X-direction)
% Y              [N_Y, N_X]   : Y matrix of the meshgrid
% O              [2 x N_scat] : Vector of the centers of the scatterers
% a              [1 x N_scat] : Vector of the radius of the scatterers
% k              [1 x 1]      : Wavenumber in the vacuum
% M_modes        [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                               number of modes kept for the Fourier basis of the
%                               scatterer number "p"
% Density        [N x 1]      : Coefficients in the Fourier basis of the
%                               density where N = sum(2*M_modes +1)
% TypeOfOperator [1 x 2]      : Coefficient of the linear combination (see
%             or [2 x N_scat]   below)
%
% Output Arguments:
% -----------------
% U              [N_Y, N_X]   : Potential computed on the grid, OUTSIDE the
%                               obstacles
%
% TypeOfOperator:
% ---------------
%
% TypeOfOperator is a 2 columns variable where (alpha_p and beta_p being 
% the scalar weight to the potentials (see above)):
%   TypeOfOperator(p, 1) = alpha_p
%   TypeOfOperator(p, 2) = beta_p
% If TypeOfOperator is a line, then:
% TypeOfOperator(p, j) = TypeOfOperator(j), forall p=1,...,N_scat and j=1,2
% -------------------------------------------------------------------------
% Example:
% --------
% if: (L_1 + 2*M_1) Density_1 + (-1)L_2 Density2
% then: TypeOfOperator = [1, 2; -1, 0]
%
% OPTIONS:
% --------
% Type: help GetPotentialOptions for more options
%
% See also InternalPotential

%%
function U = ExternalPotential(X, Y, O, a, M_modes, k, Density, TypeOfOperator, varargin)

    N_scat = length(a);
    %Get options
    [ONBOUNDARY, VERBOSITY] = GetPotentialOptions(varargin{:});
    %One density per potential ?
    SEPARATE_POTENTIAL = 0;
    if(size(Density,2) == 2)
        % There is two densities that must be computed separately 
        %(the first on single layer and the second on double layer)
       SEPARATE_POTENTIAL =1; 
    end
    
    U = zeros(size(X));

    ObstaclesMatrix = MaskMatrixObstacles(X, Y, O, a);
    if(ONBOUNDARY)
        BoundaryMatrix = BoundaryOfObstacles(X, Y, O, a);
    end
    % Initialization of the Matrix A
    % A is the matrix of the combination of the single- and double-layer 
    % integral operator (volumic)
    % The evaluation is then obtained by computing
    % U = A * rho
    Sum_modes = 2.*M_modes + 1;
    sizeA = sum(Sum_modes);
    N_x = size(X,2);
    
    %advancement counter (display purpose)
    perc_adv = 10; %percentage of advancement
    step = N_x/perc_adv; %step to reach
    adv = 1; %counter
    
    %It is possible to build one (giant) vector only, this could however
    %beat the memory limitation... A loop on every abscissa is hence prefered.
    
    for cpt_x = 1:N_x
        %extract the vector without obstacles in the column cpt_x
        GridPoints = size(ObstaclesMatrix,1)*(cpt_x-1) + find(ObstaclesMatrix(:,cpt_x)==0);
        if(ONBOUNDARY)
            GridPoints = [GridPoints; size(BoundaryMatrix,1)*(cpt_x-1) + find(BoundaryMatrix(:,cpt_x)>0)];
        end
        vect_X = X(GridPoints);
        vect_Y = Y(GridPoints);
        XY = [vect_X,vect_Y];
        if(isempty(XY))
            continue;
        end
        N_y = size(XY,1);
        %reset matrix A
        if(SEPARATE_POTENTIAL)
            ASL = zeros(N_y, sizeA); %Single layer contrib.
            ADL = zeros(N_y, sizeA); %Double layer contrib.
        else
            A = zeros(N_y, sizeA); %One matrix for the two (one density)
        end
        %Sp is a row-counter
        Sp = 0;
        %Loop on the obstacles (blocks of A)
        for p = 1:N_scat
            %Center and radius of the scatterer p
            Op = O(:,p);
            ap = a(p);
            %Vector of the modes associated with the scattererd number p
            Np = M_modes(p);
            MNp = [-Np:Np];
            matTypeOfOperator = getTypeOfOperatorAndWeight(TypeOfOperator, N_scat, p);
            if(SEPARATE_POTENTIAL)
                ASL(:,Sp + MNp +(Np+1)) = BlockPotential(XY, Op, ap, Np, k, [matTypeOfOperator(1),0], 'Outside');
                ADL(:,Sp + MNp +(Np+1)) = BlockPotential(XY, Op, ap, Np, k, [0,matTypeOfOperator(2)], 'Outside');
            else
                A(:,Sp + MNp +(Np+1)) = BlockPotential(XY, Op, ap, Np, k, matTypeOfOperator, 'Outside');
            end
            %Upgrade of the row-counter
            Sp = Sp + 2*Np+1; 
        end
        %Compute the potential on the array:
        if(SEPARATE_POTENTIAL)
            U(GridPoints) = ASL * Density(:,1) + ADL * Density(:,2);
        else
            U(GridPoints) = A * Density;
        end
        %Display advancement
        if(cpt_x > step*adv && VERBOSITY > 0)
            disp(['External Potential: ', num2str(adv*perc_adv), '% done']);
            adv = adv+1;
        end
    end
    if(VERBOSITY > 0)
        disp('External Potential: 100% done');
    end
end
%% Get the TypeOfOperator of integral operators needed
% Check if TypeOfOperator is a line or a N_scatx2 array and returns either
% the line of the line numbered p.

function [matTypeOfOperator] = getTypeOfOperatorAndWeight(TypeOfOperator, N_scat, p)
    
    if(size(TypeOfOperator, 2) ~=2)
        error('TypeOfOperator is of wrong size: must be 1x2 or N_scatx2.');
    elseif((isrow(TypeOfOperator) || iscolumn(TypeOfOperator)))
       matTypeOfOperator = TypeOfOperator;
    elseif(size(TypeOfOperator,3) == 1)
        if(size(TypeOfOperator,1) == N_scat && size(TypeOfOperator,2) == 2);
            matTypeOfOperator = TypeOfOperator(p,:);
        else
            error('TypeOfOperator is of wrong size: must be 1x2 or N_scatx2.');
        end
    else
        error('TypeOfOperator has too much dimensions!');
    end
end
