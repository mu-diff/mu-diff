% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the N_p parameters for the Fourier series truncation.
% --------------------------------------------------------------
%   M_modes = FourierTruncation(a, k)
%
% Input arguments:
% ----------------
%   
% a  [1 x N_scat]   : Radii of the circular obstacles
% k  [1 x 1]        : wavenumber in the vacuum 
% or [1 x N_scat]     or in the obstacles
%
% Output arguments:
% -----------------
% M_Modes [1 x N_scat] : Truncation index in the Fourier series of
%                        obstacles
%
% Options: 
% --------
% 1) M_modes = FourierTruncation(..., 'Min', MIN_VALUE)
% Set a minimal value for M_modes (default = 0)
%
% 2) M_modes = FourierTruncation(..., 'Tol', TOL)
% Set the tolerance in the formula to TOL
%

function M_modes = FourierTruncation(a, k, varargin)
    N_scat = length(a);
    nvarargin = length(varargin);
    cpt_arg = 1;
    MIN_VALUE = 0;
    tolM = 10.^(-10);
    
    while(cpt_arg <= nvarargin)
       if(strcmp(varargin{cpt_arg}, 'Min'))
           MIN_VALUE = varargin{cpt_arg + 1};
           cpt_arg = cpt_arg +2;
       elseif(strcmp(varargin{cpt_arg}, 'Tol'))
           tolM = varargin{cpt_arg + 1};
           cpt_arg = cpt_arg +2;
       else
          cpt_arg = cpt_arg +1; 
       end
    end
    
    M_modes = max(MIN_VALUE*ones(1,N_scat), max([2^2*ones(1,N_scat); floor(k.*a + (1/(2*sqrt(2))*log(2*sqrt(2)*pi*k.*a/tolM)).^(2/3).*(k.*a).^(1/3) +1)]));
end