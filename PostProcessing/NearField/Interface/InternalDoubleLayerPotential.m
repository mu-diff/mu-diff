% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the double-layer potential for a density "Density", on a meshgrid 
% (X, Y) and for the multiple scattering problem by circular obstacles 
% in the Fourier bases.
% -------------------------------------------------------------------------
% Note: the computation is done INSIDE the obstacles and the returned 
% variable U has 0 value OUTSIDE the obstacles.
% -------------------------------------------------------------------------
% Double-Layer potential: M lambda = sum_{p=1}^{N_scat} M_p lambda_p
% Double-Layer potential L_p rho_p associated to obstacle p, for a density lambda_p
% in H^{1/2}(Gamma_p) is given by
%    M_p lambda_p(x) = \int_{Gamma_p} dn_y G(x,y) lambda_p(y) dy, for all x in R^2\ Omega_p
%
% Omega_p: obstacle numbered p.
% Gamma_p: boundary of Omega_p.
% -------------------------------------------------------------------------
% The resulting potential is a linear combination of the single-layer 
% potentials:
%    U = sum_{p=1}^{N_scat} beta_p M_p Density_p
%
% alpha_p = scalar value
% -------------------------------------------------------------------------
% U = InternalDoubleLayerPotential(X, Y, N_scat, O, a, k, M_modes, Density, Weight)
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
% Weight         [1 x 1]      : Multiplicative coefficient (alpha_p)
%             or [N_scat x 1]   if multiple value then Weight(p) = alpha_p
%
% Output Arguments:
% -----------------
% U              [N_Y, N_X]   : Potential computed on the grid, OUTSIDE the
%                               obstacles
%
% OPTIONS:
% --------
% Type: help GetPotentialOptions for more options
%
%
% See also ExternalDoubleLayerPotential, ExternalPotential, InternalPotential

function U = InternalDoubleLayerPotential(X, Y, O, a, M_modes, k, Density, Weight, varargin)

    TypeOfOp = [0,1];
    if(isscalar(Weight))
        TypeOfOp = [Weight,0];
    elseif(iscolumn(Weight))
        TypeOfOp = [zeros(size(Weight)), Weight];
    elseif(isrow(Weight))
        TypeOfOp = [zeros(size(Weight.')), Weight.'];
    else
        error('Wrong size of Weight !');
    end
        
    U = InternalPotential(X, Y, O, a, M_modes, k, Density, TypeOfOp, varargin{:});
end