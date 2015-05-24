% mu-diff - Copyright (C) 2014-2015 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the matrix of one of the integral operator
% -------------------------------------------------------------------------
%   A = IntegralOperator(O, a, M_modes, k, TypeOfOperator)
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
%        or  [1 x N_scat]      if k is a row/column, then k(p) = wavenumber 
%                              associated to obstacle p
% TypeOfOperator (See below) : Specifies the integral operator. 
%                              See Comon/Parser.m for correspondance. 
%
% TypeOfOperator acceptable size/type:
% ------------------------------------
% - SINGLE VALUE: IntegralOperator(..., 2) or  IntegralOperator(..., {'L'})
% - ROW VECTOR/CELL: IntegralOperator(..., [1,4]) produces I+N
%                  OR IntegralOperator(..., {'I','N'}) produces I+N 
% - 2D ARRAY/CELL: IntegralOperator(..., T) then block (p,q) is of type
%   T(p,q) or T{p,q}
% - 3D ARRAY/CELL: IntegralOperator(..., T) then block (p,q) is the sum
%   of the different types: 
%       T(p,q,1) + T(p,q,2) + ... T(p,q,end)
%       T{p,q,1} + T{p,q,2} + ... T{p,q,end}
%
% OPTIONS (weight):
% -----------------
%   A = IntegralOperator(..., TypeOfOperator, Weight)
% where Weight is of the same size as TypeOfOperator and contains
% multiplicative constants to apply to the desired operator:
% - SCALAR value: IntegralOperator(..., 2, 0.5) produces 0.5*L
% - ROW VECTOR value: IntegralOperator(..., [1,4], [0.5,1]) produces 0.5*I+N
% - MATRIX value T, W: IntegralOperator(..., T, W) then block (p,q) is of type
%       W(p,q)*T(p,q)
% - 3D Matrix value T, W: IntegralOperator(..., T, W) then block (p,q) is the sum
%   of the different types: 
%       W(p,q,1)*T(p,q,1) + W(p,q,2)*T(p,q,2) + ... W(p,q,end)*T(p,q,end)
%
% EXAMPLE:
% --------
% 1) Compute N (DnSingleLayer)
%   > A = IntegralOperator(O, a, M_modes, k, 4);
%
% 2) For two obstacles in [-10,0] and [0,10] with radii 1 and 0.5, compute
% the matrix [L_11, 2*D_12; 3*N_21, 4*M_22]: 
%   > O = [-10, 0; 0, 10];
%   > a = [1, 0.5]
%   > TypeOfOp = [2,5;4,3];
%   > Weight = [1,2;3,4];
%   > A = IntegralOperator(O, a, M_modes, k, TypeOfOp, Weight);
%
% 3) Compute the single-scattering operator of the single-layer
% (off-diagonal block =0 and diagonal block = L_pp)
%   > TypeOfOp = 2*eye(N_scat, N_scat);
%   > A = IntegralOperator(O, a, M_modes, k, TypeOfOp);
%
% See also Parser, BlockIntegralOperator
% SpBlockIntegralOperator, SpIntegralOperator, SingleLayer,
% DnSingleLayer, DoubleLayer, DnDoubleLayer,
% PrecondOperator_Dirichlet, PrecondOperator_Neumann
%

function [A] = IntegralOperator(O, a, M_modes, k, TypeOfOperator, varargin)

    ScalarTypeOfOperator = ParserIntegralOperator(TypeOfOperator);
    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(size(ScalarTypeOfOperator));
    end
    
    if(size(Weight) ~=size(ScalarTypeOfOperator))
        error('Matrix of TypeOfOperator and Weight must be of same size!');
    end

    %Initialization
    N_scat = length(a);
    Sum_modes = 2.*M_modes + 1;
    sizeA = sum(Sum_modes);
    A = zeros(sizeA, sizeA);
    %Sp is a row-counter
    Sp = 0;
    %Loop on the obstacles
    for p = 1:N_scat
        %Center and radius of the scatterer p
        Op = O(:,p);
        ap = a(p);
        %Vector of the modes associated with the scattererd number p
        Np = M_modes(p); 
        MNp = [-Np:Np].';
        %Sq is a column-counter
        Sq = 0;
        %Second loop on the obstacles
        for q = 1:N_scat
            %Center and radius of the scatterer p
            Oq = O(:,q);
            aq = a(q);
            %Vector of the modes associated with the scattererd number q
            Nq = M_modes(q);
            MNq = [-Nq:Nq].';
            %Compute the block matrix
            this_TypeOfOperator = GetTypeOfOperatorOrWeight(ScalarTypeOfOperator, p, q);
            this_weight = GetTypeOfOperatorOrWeight(Weight, p, q);
            this_k = GetK(k, p);
            A(Sp + MNp +(Np+1), Sq + MNq +(Nq+1)) = ...
                BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, this_k, ...
                                      this_TypeOfOperator, this_weight, varargin{2:end});
            %Upgrade of the column-counter
            Sq = Sq + 2*Nq+1;
        end
        %Upgrade of the row-counter
        Sp = Sp + 2*Np+1; 
    end
end
