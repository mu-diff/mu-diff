% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the Far Field Pattern of a combination of single- and double-layer 
% potentials, in the case of circular scatterers.
% -------------------------------------------------------------------------
% [F] = FarField(O, a, M_modes, k, theta, Density, Weight)
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat] : Vector of the centers of the scatterers
% a       [1 x N_scat] : Vector of the radius of the scatterers
% M_modes [1 x N_Scat] : Vector of the "N_p", where "(2*N_p+1)" is the
%                        number of modes kept for the Fourier basis of the
%                        scatterer number "p".
% k            [1 x 1] : Wavenumber in the vacuum
% theta [1 x N_angles] : Angles of observation (rad).
%                        Note that theta could be a vector.
% Density      [N x 1] : Coefficients in the Fourier bases of the
%           or [N x 1]   density (or densities)
%                        N = 2*sum(M_modes+1).
% Weight       [1 x 2] : Weight to apply to the potentials (coef. alpha_p)
%      or [N_scat x 2]   if Weight is a line, then Weight(p,j) = Weight(j).
%
% Output Arguments :
% -----------------
% F     [N_angles x 1] : Vector of the Far Field F
%
% See also RCS, FarField_to_RCS
%

function F = FarField(O, a, M_modes, k, theta, Density, Weight)

    TWO_DENSITIES=0;
    if(size(Density,2)==2)
        TWO_DENSITIES =1;
    end
    N_scat = length(a);
    %thetaSER = theta*(2*pi/360);
    C = zeros(length(theta),N_scat);
    %Sp is a row-counter
    Sp = 0;
    for p=1:N_scat
        Np = M_modes(p);
        MNp = [-Np:Np];
        ap = a(p);
        Op = O(:,p);
        OO = [0;0];

        bp = norm(Op);
        alphap = fangle(Op,OO);
        EXP_1 =  exp(-1i*k*bp*cos(theta-alphap)).';
        Mde_theta1 = 1i*(theta-pi/2).'*MNp;

        matWeight = getWeight(Weight, N_scat, p);
        for iOper = 1:2
            if(matWeight(iOper) == 0)
                continue;
            end
            if(TWO_DENSITIES)
                CurrentDensity=Density(Sp + MNp +(Np+1), iOper);
            else
                CurrentDensity=Density(Sp + MNp +(Np+1));
            end
                
            this_weight = matWeight(iOper);
            FourierCoef = zeros(2*Np+1,1);
            if(this_weight ~=0)
                switch iOper
                    case 1, %Single Layer formulation
                        FourierCoef = sqrt(ap)*CurrentDensity.*besselj(MNp,k*ap).';
                    case 2,%Double Layer contribution
                        FourierCoef = -k*sqrt(ap)*CurrentDensity.*dbesselj(MNp,k*ap).';
                end
                C(:,p) = C(:,p) + this_weight*(exp(Mde_theta1)*FourierCoef).*EXP_1;
            end
            clear FourierCoef;
        end
        Sp = Sp +2*Np+1;
    end

    if (N_scat > 1)
        F = sum(1i*exp(-1i*pi/4)*sqrt(1/k)/2*C.');
    else
        F = 1i*exp(-1i*pi/4)*sqrt(1/k)/2*C.';
    end
    F = F.';
end

%% Get the Weight of integral operators needed

function [matWeight] = getWeight(Weight, N_scat, p)
    
    if(size(Weight, 2) ~=2)
        error(['Weight is of wrong size: must be 1x2 or N_scatx2.']);
    elseif((isrow(Weight) || iscolumn(Weight)))
       matWeight = Weight;
    elseif(size(Weight,3) == 1)
        if(size(Weight,1) == N_scat && size(Weight,2) == 2);
            matWeight = Weight(p,:);
        else
            error(['Weight is of wrong size: must be 1x2 or N_scatx2.']);
        end
    else
        error('Weight has too much dimensions!');
    end
end
