% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the incident wave on a grid.
% -------------------------------------------------------------------------
%   Uinc = IncidentWaveOnGrid(O, a, M_modes, k, TypeOfWave, Param)
%
% OUTPUT ARGUMETNS:
% -----------------
% Uinc [size(X)]    : Value of the incident wave on the grid
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% X           [Ny x Nx]    : Meshgrid X
% Y           [Ny x Nx]    : Meshgrid Y
% k           [1 x 1]      : Wavenumber in the vacuum
% TypeOfWave  [1 x 1]      : 'PlaneWave' (1) or 'PointSource' (2)
%            or string       
% varagin    (var)         : Parameter of the incident wave (direction, ...)
% 
% EXAMPLE:
% --------
% 1) For an incident plane wave of direction beta_inc = pi:
%   > Uinc = IncidentWaveOnGrid(X, Y, k, 1, pi)
% or
%   > Uinc = IncidentWaveOnGrid(X, Y, k, 'PlaneWave', pi)
% 2) For a point source in OS =[XS,YS]:
%   > Uinc = IncidentWaveOnGrid(X, Y, k, 2, OS)
%
% REMARK:
% -------
% - A plane wave of direction beta: uinc(X) = exp(i*beta*X)
% - A point source from OS : uinc(X) = G(X,OS) = 0.25*i*H_0^(k*\|X-OS\|)
%
% See also IncidentWave
%
%

function Uinc = IncidentWaveOnGrid(X, Y, k, TypeOfWave, varargin)
    if(isempty(varargin))
       error('Missing argument! (direction of wave or position of point source)'); 
    end
    FormattedType = ReadType(TypeOfWave);
    switch FormattedType
        case 'PlaneWave'
            beta_inc = GetParam(TypeOfWave, varargin{1});
            XYBeta = X*cos(beta_inc)+ Y*sin(beta_inc);
            Uinc = exp(1i*k*XYBeta);
        case 'PointSource'
            XS = GetParam(TypeOfWave, varargin{1});
            ABS_XY = sqrt((X-XS(1)).^2+ (Y-XS(2)).^2);
            Uinc = 1i/4*besselh(0, 1, k*ABS_XY);
        otherwise, error('Unrecognized type of incident wave');
    end
end

%%
function param = GetParam(TypeOfWave, arg)
    switch TypeOfWave
        case 'PlaneWave'
        if(isscalar(arg))
            param = arg;
        else
            error('Plane wave must be specified with a scalar argument (angle of direction)');
        end
        case 'PointSource'
            if((isrow(arg) || iscolumn(arg)) && length(arg) == 2)
                param = arg;
            else
                error('Point source must be specified with [XS,YS] argument');
            end    
    end
end

%%
function FormattedType = ReadType(TypeOfWave)
    if(ischar(TypeOfWave))
        FormattedType = TypeOfWave;
    elseif(isscalar(TypeOfWave))
        switch TypeOfWave
            case 1, FormattedType = 'PlaneWave';
            case 2, FormattedType = 'PointSource';
            otherwise, error('Unknown type of incident wave');
        end
    else
        error('Unknown type of incident wave');
    end
end