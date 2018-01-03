% mu-diff - Copyright (C) 2014-2018 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute a BLOCK matrix used to compute a linear combination of the 
% single- and double-layer potential, on a VECTOR XY and for the multiple 
% scattering problem by circular obstacles in the Fourier bases.
% -------------------------------------------------------------------------
% The computation is done INSIDE or OUTSIDE the obstacles.
% THE USER MUST HOWEVER SPECIFY IT: NO CHECK IS DONE ON THE POINTS!
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
% Ap = BlockPotential(IN_OR_OUTSIDE, XY, Op, ap, Np, kp, TypeOfOperator)
%
% Input Arguments (N_scat = number of obstacles):
% ----------------
% FORMULA        [1 x 1]      : either 'Inside' or 'Outside', this
%                               specifies the formula to use.
% XY      [N_points x 2]      : Vector of points
% Op             [2 x 1]      : Center of obstacle p
% ap             [1 x 1]      : Radius of the scatterers p
% kp             [1 x 1]      : Wavenumber to consider
% Np             [1 x 1]      : Value of "N_p", where "(2*N_p+1)" is the
%                               number of modes kept for the Fourier basis of the
%                               scatterer number "p"
% TypeOfOperator [1 x 2]      : Coefficient of the linear combination (see
%                               below)
%
% Output Arguments:
% -----------------
% Ap     [N_Points, 2*Np+1]   : Submatrix of the matrix used to compute the
%                               field on a grid
%
% TypeOfOperator:
% ---------------
%
% TypeOfOperator is a [2x1] vector where (alpha_p and beta_p being 
% the scalar weight to the potentials (see above)):
%   TypeOfOperator(1) = alpha_p
%   TypeOfOperator(2) = beta_p
%
% See also ExternalPotential, InternalPotential

%%
function Ap = BlockPotential(XY, Op, ap, Np, kp, TypeOfOperator, FORMULA)
    MNp = [-Np:Np];
    N_Points = size(XY,1);
    Ap = zeros(N_Points, 2*Np+1);

    %Rp = vector between Center Op - Points X (vector)
    % And its complex version (easier)
    Rp = XY - repeat_vert(Op.',size(XY,1));
    Rp_complex = Rp(:,1) + 1i*Rp(:,2);

    %column vector of the distance and the angle
    rp = abs(Rp_complex);
    thetap = angle(Rp_complex);

    %Trick to get a matrix with besselh and besselj and a Matlab version >= 8
    rp_mat = repeat_horiz(rp, 2*Np+1);
    MNp_mat = repeat_vert(MNp, length(rp));
    
    switch FORMULA
        case 'Inside',
            if(TypeOfOperator(1) ~= 0)
                H1m_kap = besselh(MNp_mat,1,kp*ap);
                Ap = Ap + TypeOfOperator(1)*H1m_kap;
            end
            if(TypeOfOperator(2) ~= 0)
                dH1m_kap = dbesselh(MNp_mat,1,kp*ap);
                Ap = Ap + TypeOfOperator(2)*(-kp)*dH1m_kap;
            end
            Jm_krp = besselj(MNp_mat,kp*rp_mat);
            Ap = Ap.*Jm_krp;
        case 'Outside',
            if(TypeOfOperator(1) ~= 0)
                Jm_kap = besselj(MNp_mat,kp*ap);
                Ap = Ap + TypeOfOperator(1)*Jm_kap;
            end
            if(TypeOfOperator(2) ~= 0)
                dJm_kap = dbesselj(MNp_mat,kp*ap);
                Ap = Ap + TypeOfOperator(2)*(-kp)*dJm_kap;
            end
            H1m_krp = besselh(MNp_mat,kp*rp_mat);
            Ap = Ap.*H1m_krp;
        otherwise,
            error('Unknown case');
    end
    Expthetap = exp(1i*thetap*MNp)/(sqrt(2*pi*ap));
    Ap = 1i*pi*ap/2*Ap.*Expthetap;
end
