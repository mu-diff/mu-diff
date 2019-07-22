% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Solve the Penetrable problem of multiple scattering by disk, using a 
% single-layer representation of the field outside and inside.
% -------------------------------------------------------------------------
%   Solution = PenetrableSolver(O, a, k, k_int, TypeOfWave, varargin)
%
% OUTPUT ARGUMENTS:
% -----------------
% Solution cell(7,1) : Solution such that:
%   Solution{1} = O
%   Solution{2} = a
%   Solution{3} = k
%   Solution{4} = k_int
%   Solution{5} = M_modes (Truncation index in the Fourier series of
%                              obstacles)
%   Solution{6} = TypeOfWave
%   Solution{7} = ParamOfWave (beta_inc or OS)
%   Solution{8} = rho+, the solution for the outer field
%   Solution{9} = rho-, the solution for the inner field
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
% 1- Solution = DirichletSolver(..., 'Iterative')
%  Use GMRES instead of a direct solver
% 2- Solution = DirichletSolver(..., 'Tol', TOL)
%  Set the tolerance of the iterative solver to TOL (if 'Iterative' has been
%  set)
%
% See also DirichletRCS, DirichletFarField, DirichletNearField
%

function Solution = PenetrableSolver(O, a, k, k_int, TypeOfWave, varargin)

if(isempty(varargin))
   error('Missing argument! (direction of wave or position of point source)'); 
end

%% Reading arguments
TypeOfStorage = 'Full';
TypeOfSolver = 'Direct';

nvar = length(varargin);
cpt = 2;
while (cpt < nvar)
    if(ischar(varargin{cpt}))
        if(strcmp(varargin{cpt}, 'Iterative'))
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

%% Number of obstacles
N_scat = length(a);
%% Fourier series: truncation indices
M_modes_plus = FourierTruncation(a, k, 'Min', 4);
M_modes_minus = FourierTruncation(a, k_int, 'Min', 4);
% Truncation indices are chosen as the max between intern and extern values
M_modes = max(M_modes_plus, M_modes_minus);
sum_modes = sum(2*M_modes+1);
%% GMRES parameters
RESTART = [];
TOL = 10^(-10);
MAXIT = 500;

%% -----------------------------------------
%            Integral formulation
% Single-Layer potential (plus and minus sign refer to ouside and
% inside the obstacles)
% u = \int_{\Gamma} G^{+/-}(x,y) rho^{+/-}(y) d\Gamma(y),  x
% rho^{+/-} = densities (unknown)
% vector of the "lambda" (contrast)
%% -----------------------------------------
%% Boundary integral operators
%Single-layer (exterior)
Splus = SingleLayer(O, a, M_modes, k); 
%Single-layer (intern: no off-diagonal blocks but different wave numbers)
Sminus = IntegralOperator(O, a, M_modes, k_int, 2*eye(N_scat,N_scat));
%Dn Single Layer (exterior)
Nplus = DnSingleLayer(O, a, M_modes, k);
%Dn Single Layer (intern, no off-diagonal blocks but different wave numbers)
Nminus = IntegralOperator(O, a, M_modes, k_int, 4*eye(N_scat,N_scat));
%identity matrix
Identity = eye(size(Nplus));
%% Boundary integral operator
A = [Splus, -Sminus; -0.5*Identity + Nplus, -(0.5*Identity + Nminus)];
%cleaning memory
clear Splus Sminus Nplus Nminus Identity;
%% Incident plane wave (right hand side)
ParamWave = varargin{1};
Uinc = IncidentWave(O, a, M_modes, k, TypeOfWave, ParamWave);
DnUinc = IncidentWave(O, a, M_modes, k, ['Dn',TypeOfWave], ParamWave);
[B] = [Uinc;DnUinc];
clear Uinc DnUinc;
%% Solving linear system
if(strcmp(TypeOfSolver, 'Direct'))
    rho = A\B;
else
    rho = gmres(A, B, RESTART, TOL, MAXIT, [], []);
end        
%% Extracting the "exterior" and "interior" densities
rho_plus = rho(1:sum_modes);
rho_minus = rho(sum_modes+1:end);
clear rho;
%% SOLUTION
Solution = cell(7,1);

Solution{1} = O;
Solution{2} = a;
Solution{3} = k;
Solution{4} = k_int;
Solution{5} = M_modes;
Solution{6} = TypeOfWave;
Solution{7} = ParamWave;
Solution{8} = rho_plus;
Solution{9} = rho_minus;


