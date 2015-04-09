% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the block vector of (the opposite of) the coefficients of an 
% incident wave, either the trace of the normal derivative trace, on one
% of the obstacles, in the Fourier bases.
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -u^inc, -dn
% u^inc, ...
% -------------------------------------------------------------------------
%   Bp = BlockIncidentWave(Op, ap, Np, k, TypeOfWave, Param)
%
% OUTPUT ARGUMETNS:
% -----------------
% B [2*Np+1,1]    : Vector of the coefficients
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% Op          [2 x N_scat] : Vector of the center of the scatterer
% ap          [1 x N_scat] : Radius of the scatterer
% Np          [1 x N_scat] : Truncation index in the Fourier series of
%                            obstacles
% k           [1 x 1]      : Wavenumber in the vacuum
% TypeOfWave  [1 x 1]      : Plane wave (1), Dn Plane wave (2), 
%                            Point source (3), Dn Point source (4),
%                            Precond Plane wave (5), Precond Dn Plane wave (6)
% Param      (var)         : Parameter of the incident wave (direction, ...)
% 
% REMARK:
% -------
% - It is not possible to mix plane wave and point source. This should be
% done by two call of this functions.
% - A plane wave of direction beta: uinc(X) = exp(i*beta*X)
% - A point source from OS : uinc(X) = G(X,OS) = 0.25*i*H_0^(k*\|X-OS\|)
%
% See also Incidentwave, PlaneWave, DnPlaneWave, PointSource, DnPointSource,
% PlaneWavePrecond, DnPlaneWavePrecond
%
%

function Bp = BlockIncidentWave(Op, ap, Np, k, TypeOfWave, varargin)

    if(isempty(varargin))
       error('Missing argument! (direction of wave or position of point source)'); 
    end
    if(~isscalar(TypeOfWave))
       error('In block function, TypeOfWave must be a scalar');
    end
    Param = varargin{1};
    %Initialization
    Bp = zeros(2*Np+1,1);
    MNp = [-Np:Np].';
    ScalarTypeOfWave = ParserIncidentWave(TypeOfWave);
    switch ScalarTypeOfWave
        case 0, %Nothing
        case 1, %Plane wave
            beta_inc = GetParam(TypeOfWave, Param);
            Jm_kap = besselj(MNp,k*ap);
            I = exp(1i*k*(Op(1)*cos(beta_inc) + Op(2)*sin(beta_inc)));
            d_m = I*exp(1i*MNp*(pi/2-beta_inc));
            Bp = - sqrt(2*pi*ap)* (d_m.*Jm_kap);                
        case 2, %Dn Plane wave
            beta_inc = GetParam(TypeOfWave, Param);
            dJm_kap = dbesselj(MNp,k*ap);
            I = exp(1i*k*(Op(1)*cos(beta_inc) + Op(2)*sin(beta_inc)));
            d_m = I*exp(1i*MNp*(pi/2-beta_inc));
            Bp= - k*sqrt(2*pi*ap)* (d_m.*dJm_kap);
        case 3, %Point source
            XS = GetParam(TypeOfWave, Param);
            betap = fangle(Op,XS);
            bp = norm(Op-XS);
            Jm_kap = besselj(MNp,k*ap);
            Hm_kbp = besselh(-MNp,1,k*bp);
            Exp_betap = exp(-1i*MNp*betap);    
            Bp= - 1i/4*sqrt(2*pi*ap)*Jm_kap.*Hm_kbp.*Exp_betap;
        case 4, %Dn Point source
            XS = GetParam(TypeOfWave, Param);
            betap = fangle(Op,XS);
            bp = norm(Op-XS);
            dJm_kap = dbesselj(MNp,k*ap);
            Hm_kbp = besselh(-MNp,1,k*bp);
            Exp_betap = exp(-1i*MNp*betap);
            Bp= - 1i*k/4*sqrt(2*pi*ap)*dJm_kap.*Hm_kbp.*Exp_betap;
        case 5, %Precond Dirichlet Plane Wave (trace precond by single-layer)
            beta_inc = GetParam(TypeOfWave, Param);
            H1m_kap_inv = 1./besselh(MNp,1,k*ap);
            I = exp(1i*k*(Op(1)*cos(beta_inc) + Op(2)*sin(beta_inc)));
            d_m = I*exp(1i*MNp*(pi/2-beta_inc));
            Bp= 1i*2*sqrt(2/pi/ap)* (d_m.*H1m_kap_inv);
        case 6, %Precond Neumann Dn Plane Wave (dn trace precond by dn double-layer)
            beta_inc = GetParam(TypeOfWave, Param);
            dH1m_kap_inv = 1./dbesselh(MNp,1,k*ap);
            I = exp(1i*k*(Op(1)*cos(beta_inc) + Op(2)*sin(beta_inc)));
            d_m = I*exp(1i*MNp*(pi/2-beta_inc));
            Bp= 1i*2*sqrt(2/pi/ap)/k* (d_m.*dH1m_kap_inv);                
        case 7, %Precond Dirichlet Point Source (trace precond by single-layer)
            XS = GetParam(TypeOfWave, Param);
            betap = fangle(Op,XS);
            bp = norm(Op-XS);
            Hm_kbp = besselh(-MNp,1,k*bp);
            H1m_kap_inv = 1./besselh(MNp,1,k*ap);
            Exp_betap = exp(-1i*MNp*betap);    
            Bp= - 1/sqrt(2*pi*ap)*Hm_kbp.*H1m_kap_inv.*Exp_betap;
        case 8, %Precond Neumann Dn Point Source (dn trace precond by double-layer)
            XS = GetParam(TypeOfWave, Param);
            betap = fangle(Op,XS);
            bp = norm(Op-XS);
            Hm_kbp = besselh(-MNp,1,k*bp);
            dH1m_kap_inv = 1./dbesselh(MNp,1,k*ap);
            Exp_betap = exp(-1i*MNp*betap);
            Bp= 1/(k*sqrt(2*pi*ap))*Hm_kbp.*dH1m_kap_inv.*Exp_betap;
    end
end

%%
function param = GetParam(TypeOfWave, arg)
    if( TypeOfWave == 1 || TypeOfWave == 2 || (TypeOfWave >=5 && TypeOfWave <=6))
        %Plane wave
        if(isscalar(arg))
            param = arg;
        else
            error('Plane wave must be specified with a scalar argument (angle of direction)');
        end
    else
         %Point source
        if((isrow(arg) || iscolumn(arg)) && length(arg) == 2)
            if(isrow(arg))
                param = arg.';
            else
                param = arg;
            end
        else
            error('Point source must be specified with [XS,YS] argument');
        end    
    end
end