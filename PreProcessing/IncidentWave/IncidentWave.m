% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the vector of (the opposite of) the coefficients of an 
% incident wave, either the trace of the normal derivative trace, for the 
% multiple scattering problem by disks, in the Fourier bases.
% -------------------------------------------------------------------------
% REMARK: What is computed is the OPPOSITE of the coefficient: -u^inc, -dn
% u^inc, ...
% -------------------------------------------------------------------------
%   B = IncidentWave(O, a, M_modes, k, TypeOfWave, Param)
%
% OUTPUT ARGUMETNS:
% -----------------
% B [sum(2*M_modes+1),1]    : Vector of the coefficients
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% O           [2 x N_scat] : Vector of the centers of the scatterers
% a           [1 x N_scat] : Vector of the radius of the scatterers
% M_Modes     [1 x N_scat] : Truncation index in the Fourier series of
%                            obstacles
% k           [1 x 1]      : Wavenumber in the vacuum
% TypeOfWave  [1 x 1]      : Plane wave (1), Dn Plane wave (2), 
%          or [N_scat x 1]   Point source (3), Dn Point source (4),
%                            Precond Plane wave (5), Precond Dn Plane wave (6)
% varagin    (var)         : Parameter of the incident wave (direction, ...)
% 
% EXAMPLE:
% --------
% 1) For an incident plane wave of direction beta_inc = pi:
%   > B = IncidentWave(O, a, M_modes, k, 1, beta_inc)
% 2) For a point source in OS =[XS,YS]:
%   > B = IncidentWave(O, a, M_modes, k, 3, OS)
% 3) For 2 obstacles, a vector with (-u^inc, -dn u^inc) where u^inc is a
% plane wave of direction beta_inc:
%   > B = IncidentWave(O, a, M_modes, k, [1,2], beta_inc)
%
% REMARK:
% -------
% - It is not possible to mix plane wave and point source. This should be
% done by two call of this functions.
% - A plane wave of direction beta: uinc(X) = exp(i*beta*X)
% - A point source from OS : uinc(X) = G(X,OS) = 0.25*i*H_0^(k*\|X-OS\|)
%
% See also PlaneWave, DnPlaneWave, PointSource, DnPointSource,
% PlaneWavePrecond, DnPlaneWavePrecond
%
%
%%
function B = IncidentWave(O, a, M_modes, k, TypeOfWave, varargin)

    if(isempty(varargin))
       error('Missing argument! (direction of wave or position of point source)'); 
    end
    %Initialization
    N_scat = length(a);
    Sum_modes = 2.*M_modes + 1;
    B = zeros(sum(Sum_modes),1);
    %Sp is a row-counter
    Sp = 0;
    %Loop on the obstacles (blocks of u_inc)
    for p = 1:N_scat
        %Center and radius of the scatterer p
        Op = O(:,p);
        ap = a(p);
        %Vector of the modes associated with the scattererd number p
        Np = M_modes(p);
        MNp = [-Np:Np].';
        
        this_wave = GetTypeOfWave(TypeOfWave, p);
        B(Sp + MNp +(Np+1))= BlockIncidentWave(Op, ap, Np, k, this_wave, varargin{1});
        %Upgrade of the row-counter
        Sp = Sp + 2*Np+1; 
    end
end

%%
function this_wave = GetTypeOfWave(TypeOfWave, p)
    if(length(TypeOfWave) == 1)
        this_wave= TypeOfWave;
    else
        this_wave = TypeOfWave(p);
    end

end

%%
function param = GetParam(this_wave, arg)
    if( this_wave == 1 || this_wave == 2 || (this_wave >=5 && this_wave <=6))
        %Plane wave
        if(isscalar(arg))
            param = arg;
        else
            error('Plane wave must be specified with a scalar argument (angle of direction)');
        end
    else
         %Point source
        if((isrow(arg) || iscolumn(arg)) && length(arg) == 2)
            param = arg;
        else
            error('Point source must be specified with [XS,YS] argument');
        end    
    end
end           
