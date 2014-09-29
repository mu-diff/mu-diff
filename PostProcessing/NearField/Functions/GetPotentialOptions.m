% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Get the options for the potentials functions on near field computation
% Every function involving the computation of a near field potential have
% the following option ("Fun" being a generic function):
%
%   Fun(..., 'Verbosity', VERBOSITY)
% VERBOSITY [1x1] : set to 0 to not avoid message
%
%   Fun(..., 'OnBoundary', ONBOUNDARY)
% Compute the solution also on the boundaries of the obstacles
% (ONBOUNDARY = 1) or not (ONBOUNDARY = 0, default value)
%
%
% See also ExternalPotential, InternalPotential
%
function [ONBOUNDARY, VERBOSITY] = GetPotentialOptions(varargin)

    ONBOUNDARY = 0;
    VERBOSITY = 1;

    nvarargin = length(varargin);
    cpt_arg = 1;
    while(cpt_arg <= nvarargin)
           if(strcmp(varargin{cpt_arg}, 'Verbosity'))
               if(isscalar(varargin{cpt_arg+1}))
                   VERBOSITY = varargin{cpt_arg+1};
                   cpt_arg = cpt_arg +2;
               else
                   error('Verbosity must be scalar');
               end
           elseif(strcmp(varargin{cpt_arg}, 'OnBoundary'))           
               ONBOUNDARY = varargin{cpt_arg+1};
               cpt_arg = cpt_arg + 2;
           else % next one
               cpt_arg = cpt_arg +1;
           end
    end
end