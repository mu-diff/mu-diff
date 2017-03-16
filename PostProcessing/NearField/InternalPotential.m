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
% The computation is done INSIDE the obstacles and the returned variable U
% outside the obstacles is set to 0.
% -------------------------------------------------------------------------
% Single-Layer potential L_p rho_p associated to obstacle p, for a density rho_p
% in H^{-1/2}(Gamma_p) is given by
%    L_p rho_p(x) = \int_{Gamma_p} G_p(x,y) rho_p(y) dy, for all x in Omega_p
% Double-Layer potential M_p lambda_p for a density lambda_p in H^{1/2}(Gamma_p):
%    M lambda(x) = - \int_{Gamma_p} dn_y G_p(x,y) lambda_p(y) dy, for all x in Omega_p
%
% Omega_p: obstacle numbered p.
% Gamma_p: boundary of Omega_p.
% Note: the wavenumber contained in the Green function G_p may changed from
%       one obstacle to another
% -------------------------------------------------------------------------
% The resulting potential is a linear combination of the single- and
% double-layer potentials, but associated to the SAME potential Density_p:
%    U = sum_{p=1}^{N_scat} (alpha_p L_p + beta_p M_p) Density_p
% -------------------------------------------------------------------------
% U = InternalPotential(X, Y, N_scat, O, a, k, M_modes, Density, TypeOfOperator)
%
% Input Arguments (N_scat = number of obstacles):
% ----------------
% X              [N_Y, N_X]   : X matrix of the meshgrid (N_Y = Nb. points
%                               in Y-direction and N_X = Nb. points in X-direction)
% Y              [N_Y, N_X]   : Y matrix of the meshgrid
% O              [2 x N_scat] : Vector of the centers of the scatterers
% a              [1 x N_scat] : Vector of the radius of the scatterers
% k_int          [1 x 1]      : Wavenumber in the obstacle
%             or [1 x N_scat]   if k_int is not a scalar then k inside obstacle p
%                               is given by k_int(p)
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
% U              [N_Y, N_X]   : Potential computed on the grid, INDIDE the
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
% if: (L_1 + 2*M_1) Density_1 + (-1)L_2 Density2
% then: TypeOfOperator = [1, 2; -1, 0]
%
% OPTIONS:
% --------
% Type: help GetPotentialOptions for more options
%
% See also ExternalPotential

%%
function U = InternalPotential(X, Y, O, a, M_modes, k_int, Density, TypeOfOperator, varargin)


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
    %Sum of the modes
    Sum_modes = 2.*M_modes + 1;
    sizeA = sum(Sum_modes);
    %counter: where to begin in vector Density ?
    Sp = 0;
    for p=1:N_scat % obstacle where to compute 
        %Building one vector to reduce the computation to one matrix-vector
        %product thanks to a the mask matrix
        ObstacleMatrix_p = find(ObstaclesMatrix == p);
        if(ONBOUNDARY)
            ObstacleMatrix_p = [ObstacleMatrix_p; find(ObstaclesMatrix == p+0.5)];
        end
        vect_Xq = X(ObstacleMatrix_p);
        vect_Yq = Y(ObstacleMatrix_p);
        XYp = [vect_Xq,vect_Yq];
        if(isempty(XYp))
            continue;
        end
        %reset matrix A
        if(SEPARATE_POTENTIAL)
            ASL = zeros(size(XYp,1), sizeA); %Single Layer contrib.
            ADL = zeros(size(XYp,1), sizeA); %Single Layer contrib.
        else
            A = zeros(size(XYp,1), sizeA); %Single- and Double- Layer contrib.
        end
        %Center and radius of the scatterer p
        Op = O(:,p);
        ap = a(p);
        %Vector of the modes associated with the scatterer numbered p
        Np = M_modes(p);
        MNp = [-Np:Np];
        %Interne wavenumber
        if(isscalar(k_int))
            kp = k_int;
        else
            kp = k_int(p);
        end
        matTypeOfOperator = getTypeOfOperatorAndWeight(TypeOfOperator, N_scat, p);
        if(SEPARATE_POTENTIAL)
            ASL = BlockPotential(XYp, Op, ap, Np, kp, [matTypeOfOperator(1), 0], 'Inside');
            ADL = BlockPotential(XYp, Op, ap, Np, kp, [0, matTypeOfOperator(2)], 'Inside');
        else
            A = BlockPotential(XYp, Op, ap, Np, kp, matTypeOfOperator, 'Inside');
        end
        %Compute the potential inside obstacle p
        if(SEPARATE_POTENTIAL)
            U(ObstacleMatrix_p) = ASL*Density(Sp + MNp +(Np+1),1) + ADL*Density(Sp + MNp +(Np+1),2);
        else
            U(ObstacleMatrix_p) = A*Density(Sp + MNp +(Np+1),1);
        end
        %Display advancement
        if(VERBOSITY > 0)
            disp(['Internal Potential, obstacle ', num2str(p), ': done']);
        end
        %Upgrade of the row-counter
        Sp = Sp + 2*Np+1;
    end
end
%% Get the TypeOfOperator of integral operators needed
% Check if TypeOfOperator is a line or a N_scatx2 array and returns either
% the line of the line numbered p.

function [matTypeOfOperator] = getTypeOfOperatorAndWeight(TypeOfOperator, N_scat, p)
    
    if(size(TypeOfOperator, 2) ~=2)
        error(['TypeOfOperator is of wrong size: must be 1x2 or N_scatx2.']);
    elseif((isrow(TypeOfOperator) || iscolumn(TypeOfOperator)))
       matTypeOfOperator = TypeOfOperator;
    elseif(size(TypeOfOperator,3) == 1)
        if(size(TypeOfOperator,1) == N_scat && size(TypeOfOperator,2) == 2);
            matTypeOfOperator = TypeOfOperator(p,:);
        else
            error(['TypeOfOperator is of wrong size: must be 1x2 or N_scatx2.']);
        end
    else
        error('TypeOfOperator has too much dimensions!');
    end
end
