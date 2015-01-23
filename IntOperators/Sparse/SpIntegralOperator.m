% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the matrix of one of the integral operator
% SPARSE VERSION
% -------------------------------------------------------------------------
%   A = SpIntegralOperator(O, a, M_modes, k, TypeOfOperator)
%
% OUTPUT ARGUMENTS:
% -----------------
% A [2*sum(M_modes+1) x 2*sum(M_modes+1)] : Matrix of the system
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------
% O          [2 x N_scat]    : Coordinates of the center of the disk 1
% a          [1 x N_scat]    : Radius of disk 1
% M_Modes    [1 x N_scat]    : Truncation index in the Fourier series of
%                              obstacles
% k          [1 x 1]         : Wavenumber in the vacuum
% TypeOfOperator  (see below): Null matrix (0 or 'Z'), Identity I (1 or 'I'), 
%                              SingleLayer L (2 or 'L'), DoubleLayer M (3 or 'M'),
%                              DnSingleLayer N (4 or 'N'), DnDoubleLayer D (5 or 'D')
%                              Precond_Dirichlet (6 or 'P'), Precond_Neumann (7 or 'Q')
% 
% TypeOfOperator:
% ---------------
% - integer or char values
% - SCALAR value: SpIntegralOperator(..., 2) produces SingleLayer L
%      OR one char value: SpIntegralOperator(..., 'L') produces SingleLayer L
% - MATRIX value T: IntegralOperator(..., T) then block (p,q) is of type
%   T(p,q)
%
% OPTIONS (weight):
% -----------------
%   A = IntegralOperator(..., TypeOfOperator, Weight)
% where Weight is of the same size as TypeOfOperator and contains
% multiplicative constants to apply to the desired operator:
% - SCALAR value: SpIntegralOperator(..., 2, 0.5) produces 0.5*L
% - MATRIX value T, W: SpIntegralOperator(..., T, W) then block (p,q) is of type
%       W(p,q)*T(p,q)
%
% REMARK: contrary to the dense function, it is not possible to directly sum
% two operators in sparse version. This can be done however in the Matrix-vector
% product functon.
%
% EXAMPLE:
% --------
% 1) Compute N (DnSingleLayer)
%   > A = cell(3,1);
%   > [A{1},A{2},A{3}] = SpIntegralOperator(O, a, M_modes, k, 4)
% 2) For two obstacles in [-10,0] and [0,10] with radii 1 and 0.5, compute
% the matrix [L_11, 2*D_12; 3*N_21, 4*M_22]: 
%   > O = [-10, 0; 0, 10];
%   > a = [1, 0.5]
%   > A = cell(3,1);
%   > TypeOfOp = [2,5;4,3];
%   > Weight = [1,2;3,4];
%   > [A{1},A{2},A{3}] = SpIntegralOperator(O, a, M_modes, k, TypeOfOp, Weight)
%
% See also Parser, BlockIntegralOperator
% SpBlockIntegralOperator, SpIntegralOperator, SpSingleLayer,
% SpDnSingleLayer, SpDoubleLayer, SpDnDoubleLayer,
% SpPrecondOperator_Dirichlet, SpPrecondOperator_Neumann

%%
function A = SpIntegralOperator(O, a, M_modes, k, TypeOfOperator, varargin)
    TypeOfOperator = IntegralOperatorParser(TypeOfOperator);
    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(size(TypeOfOperator));
    end
    
    if(size(Weight) ~=size(TypeOfOperator))
        error('Matrix of TypeOfOperator and Weight must be of same size!');
    end

    N_scat = length(a);
    %initialization of the Matrix A
    Nmax = max(M_modes);
    SizeMax = 2*Nmax+1;
    A = cell(3,1);
    A{1} = zeros(SizeMax, N_scat, N_scat);
    A{2} = zeros(2*SizeMax-1, N_scat, N_scat);
    A{3} = zeros(SizeMax, N_scat, N_scat);

    %Loop on the obstacles (blocks of A)
    for p = 1:N_scat
        %Center and radius of the scatterer p
        Op = O(:,p);
        ap = a(p);
        Np = M_modes(p);
        for q = 1:N_scat
            %Center and radius of the scatterer p
            Oq = O(:,q);
            aq = a(q);
            Nq = M_modes(q);
            %Compute the block matrix
            [this_type, this_weight] = getTypeOfOperatorAndWeight(TypeOfOperator, Weight, p, q);
            [A{1}(:,p,q), A{2}(:,p,q), A{3}(:,p,q)] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, k, this_type, this_weight);
        end
    end
end

%% Get the TypeOfOperator of integral operators needed and (if any) the associated weight

function [this_type, this_weight] = getTypeOfOperatorAndWeight(TypeOfOperator, Weight, p, q)
    
    if (isscalar(TypeOfOperator))
       this_type = TypeOfOperator;
    elseif(size(TypeOfOperator,1) == size(TypeOfOperator,2))
       this_type = TypeOfOperator(p,q);
    else
        error('Wrong size of TypeOfOperator (either a scalar or a N_scatxN_scat matrix)');
    end
    
    if (isscalar(Weight))
       this_weight = Weight;
    elseif(size(Weight,1) == size(Weight,2))
       this_weight = Weight(p,q);
    else
        error('Wrong size of TypeOfOperator (either a scalar or a N_scatxN_scat matrix)');
    end
    
end
