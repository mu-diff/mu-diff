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
% TypeOfOperator (See below) : Specifies the integral operator. 
%                              See Comon/Parser.m for correspondance. 
%
% TypeOfOperator acceptable size/type:
% ------------------------------------
% - SINGLE VALUE: SpIntegralOperator(..., 2) or   SpIntegralOperator(..., {'L'})
% - 2D ARRAY/CELL: IntegralOperator(..., T) then block (p,q) is of type
%   T(p,q) or T{p,q}
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
    N_scat = length(a);
    CheckSize(N_scat, TypeOfOperator, 'TypeOfOperator');    

    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(size(TypeOfOperator));
    end
    CheckSize(N_scat, Weight, 'Weight');
    if(size(Weight) ~=size(TypeOfOperator))
        error('Matrix of TypeOfOperator and Weight must be of same size!');
    end

    %Parse the typeofop to get a purely numeric array
    TypeOfOperator = Parser(TypeOfOperator);

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
            %Get the type of the block and the weight
            this_type = GetTypeOfOperatorOrWeight(TypeOfOperator, p, q);
            this_weight = GetTypeOfOperatorOrWeight(Weight, p, q);
            this_k = GetK(k, p);
            %Compute the block matrix
            [A{1}(:,p,q), A{2}(:,p,q), A{3}(:,p,q)] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, this_k, this_type, this_weight);
        end
    end
end

%%
function []=  CheckSize(N_scat, A, name)
    ni = size(A,1);
    nj = size(A,2);
    nk = size(A,3);
    isErr = false;
    if(nk>1)
        isErr = true;
    end
    if(ni ~= 1 && ni ~= N_scat)
        isErr = true;
    end
    if(ni ~=nj)
        isErr = true;
    end
    if(isErr)
            error('Wrong size of ', name, [' (either one valued or a ' ...
                            'N_scat x N_scat array or cell)']);
    end
end