% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Solve the Neumann problem of multiple scattering by disk, using a 
% double-layer representation of the scattered field.
% -------------------------------------------------------------------------
%   Solution = NeumannSolver(O, a, k, TypeOfWave, ParamWave)
%
% OUTPUT ARGUMENTS:
% -----------------
% Solution cell(7,1) : Solution such that:
%   Solution{1} = O
%   Solution{2} = a
%   Solution{3} = k
%   Solution{4} = M_modes (Truncation index in the Fourier series of
%                              obstacles)
%   Solution{5} = TypeOfWave
%   Solution{6} = ParamOfWave (beta_inc or OS)
%   Solution{7} = lambda, the solution
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% O          [2 x N_scat]    : Coordinates of the center of the disk 1
% a          [1 x N_scat]    : Radius of disk 1
% k          [1 x 1]         : Wavenumber in the vacuum
% TypeOfWave [char]          : 'PlaneWave' or 'PointSource'
% ParamWave  [1 x 1]         : either the angle of incidence (PlaneWave)
%              or [2 x 1]      or the coordinates of the point source
%
% OPTIONS :
% ----------
% 1- Solution = NeumannSolver(..., 'Sparse')
%  Force sparse storage for solving process (otherwise, Full storage)
%  This should be used if the problem is large.
% 2- Solution = NeumannSolver(..., 'Iterative')
%  Use GMRES instead of a direct solver
%  This should be used if the problem is large.
% 3- Solution = NeumannSolver(..., 'Tol', TOL)
%  Set the tolerance of the iterative solver to TOL (if 'Iterative' has been
%  set). Default: 1e-6
%
% See also NeumannRCS, NeumannFarField, NeumannNearField
%
function Solution = NeumannSolver(O, a, k, TypeOfWave, varargin)
if(isempty(varargin))
   error('Missing argument! (direction of wave or position of point source)'); 
end

%% Fourier series: truncation indices
M_modes = FourierTruncation(a, k, 'Min', 1);
%% GMRES parameters
MAXIT = sum(2*M_modes+1);
RESTART = [];
TOL = 10^(-6);
%% Reading arguments
TypeOfStorage = 'Full';
TypeOfSolver = 'Direct';

nvar = length(varargin);
cpt = 2;
while (cpt < nvar)
    if(ischar(varargin{cpt}))
        if(strcmp(varargin{cpt}, 'Sparse'))
            TypeOfStorage = 'Sparse';
            cpt = cpt + 1;
        elseif(strcmp(varargin{cpt}, 'Iterative'))
            TypeOfSolver =  'Iterative';
            cpt = cpt + 1;
        elseif(strcmp(varargin{cpt}, 'Tol'))
            TOL =  varargin{cpt+1};
            cpt = cpt + 2;
        else
            warning('Unknown option')
            cpt = cpt + 1;
        end
    else
        warning('Unknown option')
        cpt = cpt + 1;
    end
end

if(strcmp(TypeOfStorage, 'Sparse'))
    TypeOfSolver = 'Iterative';
    disp('Warning: Direct solver cannot be used with sparse storage');
end

%% Integral formulations
%% Single scattering preconditioned integral equation (double-layer formulation)
% Incident plane wave (right hand side)
ParamWave = varargin{1};
TypeOfWavePrecond = ['Dn', TypeOfWave, 'Precond'];
DnUincPrecond = IncidentWave(O, a, M_modes, k, TypeOfWavePrecond, ParamWave);
%MATRIX and SOLUTION
switch TypeOfStorage
    case 'Full',
        A_Precond = IntegralOperator(O, a, M_modes, k, {'Dprec'});
        if(strcmp(TypeOfSolver, 'Direct'))
            lambda = A_Precond\DnUincPrecond;
        else
            lambda = gmres(A_Precond, DnUincPrecond, RESTART, TOL, MAXIT, [], []);
        end        
    case 'Sparse',
        SpPrecond = SpIntegralOperator(O, a, M_modes, k, {'Dprec'});
        lambda = gmres(@(X)SpMatVec(X,M_modes,SpPrecond), DnUincPrecond, RESTART, TOL, MAXIT, [], []);
end

%% SOLUTION
Solution = cell(7,1);
Solution{1} = O;
Solution{2} = a;
Solution{3} = k;
Solution{4} = M_modes;
Solution{5} = TypeOfWave;
Solution{6} = ParamWave;
Solution{7} = lambda;

