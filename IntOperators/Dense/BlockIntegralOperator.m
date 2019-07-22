% mu-diff - Copyright (C) 2014-2019 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the interaction matrix of one of the integral operator in 
% the Fourier basis
% -------------------------------------------------------------------------
% Apq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, TypeOfOperator, ...)
%
% OUTPUT ARGUMENTS:
% -----------------
% Apq  [2*Np+1 x 2*Nq+1]  : submatrix of the system
%
% INPUT ARGUMENTS:
% ----------------
% Op              [2 x 1]    : Coordinates of the center of the disk p
% ap              [1 x 1]    : Radius of disk p
% Np              [1 x 1]    : Truncation index in the Fourier series of
%                              obstacle p
% Oq              [2 x 1]    : Coordinates of the center of the disk q
% aq              [1 x 1]    : Radius of disk q
% Nq              [1 x 1]    : Truncation index in the Fourier series of 
%                              obstacle q
% k               [1 x 1]    : Wavenumber in the vacuum
% TypeOfOperator (See below) : Specifies the integral operator. 
%                              See Comon/Parser.m for correspondance. 
%
% TypeOfOperator acceptable size/type:
% ------------------------------------
% - SINGLE VALUE: BlockIntegralOperator(..., 2) or BlockIntegralOperator(..., {'L'})
% - ROW VECTOR/CELL: multiple operators are computed and  together. For example:
%   TypeOfOperator = [1,4] or TypeOfOperator = {'I','N'}
%     will produce I + N
% 
% OPTIONS:
% --------
% Apq = blockIntegralOperator(..., TypeOfOperator, Weight)
% Weight  [1 x 1]    : multiplicative constant to apply to the operator(s)
%    or [1 x Nop]    
%
% EXAMPLES:
% ---------
% BlockIntegralOperator(..., 2, 0.5) will produce 0.5*L
%   OR for the same result: BlockIntegralOperator(..., 'L', 0.5)
% BlockIntegralOperator(..., [1,4], [0.5,1]) will produce: 0.5*I + N
%   OR BlockIntegralOperator(..., {'L','N'}, [0.5,1]) will produce: 0.5*I + N
% 
% See also Parser, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator, BlockSingleLayer,
% BlockDnSingleLayer, BlockDoubleLayer, BlockDnDoubleLayer,
% BlockPrecondDirichlet, BlockPrecondNeumann

%%
function Apq = BlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, k, TypeOfOperator, varargin)
    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(size(TypeOfOperator));
    end
    if(size(Weight) ~=size(TypeOfOperator))
        error('Matrix of TypeOfOperator and Weight must be of same size!');
    end

    ScalarTypeOfOperator = ParserIntegralOperator(TypeOfOperator);
    nTypeOfOperator = length(ScalarTypeOfOperator);

    %Initialization of the submatrix Apq
    Apq = zeros(2*Np+1, 2*Nq+1);

    %Vector of the modes associated with the scatterer 1 and 2
    MNp = [-Np:Np].';

    if (Op(1) == Oq(1) && Op(2) == Oq(2) && ap == aq) % Diagonal block (which is diagonal)
        if(Np ~= Nq)
            error('Obstacles are the same but with different number of modes !')
        end
        for iOperator=1:nTypeOfOperator
            this_Operator = ScalarTypeOfOperator(iOperator);
            this_weight = Weight(iOperator);
            if(this_weight == 0)
                continue;
            end
            switch this_Operator
                case -1, %Custom operator
                    Apq = Apq + this_weight*varargin{2}(Op, ap, Np, Oq, aq, Nq, k, varargin{3:end});
                case 0, %null Matrix
                case 1, % Identity
                    Apq = Apq + this_weight*eye(2*Np+1, 2*Np+1);
                case 2, % Single Layer
                    Jm_kap = besselj(MNp,k*ap);
                    H1m_kap = besselh(MNp,1,k*ap);
                    Apq = Apq + this_weight*diag(1i*ap*pi*0.5*H1m_kap.*Jm_kap);
                case 3, % Double Layer
                    dJm_kap = dbesselj(MNp,k*ap);
                    H1m_kap = besselh(MNp,1,k*ap);
                    Apq = Apq + this_weight*diag(0.5 - 1i*k*ap*pi*0.5*H1m_kap.*dJm_kap);
                case 4, % Dn Single Layer
                    dJm_kap = dbesselj(MNp,k*ap);
                    H1m_kap = besselh(MNp,1,k*ap);
                    Apq = Apq + this_weight*diag(-0.5 + 1i*k*ap*pi*0.5*H1m_kap.*dJm_kap);
                case 5, % Dn Double Layer
                    dJm_kap = dbesselj(MNp,k*ap);
                    dH1m_kap = dbesselh(MNp,1,k*ap);
                    Apq = Apq + this_weight*diag(- 1i*k^2*ap*pi*0.5*dH1m_kap.*dJm_kap);
                case 6, % Single scattering preconditioned integral operator for Dirichlet
                    Apq = Apq + this_weight*eye(2*Np+1, 2*Np+1);
                case 7, % Single scattering preconditioned integral operator for Neumann
                    Apq = Apq + this_weight*eye(2*Np+1, 2*Np+1);
            end
        end
    else %Off-diagonal blocks
        if(norm(Op-Oq) - ap -aq <= 0)
           error('Obstacles are overlapping !'); 
        end
        MNq = [-Nq:Nq].';
        % Distance and angle between centers
        bpq = norm(Op - Oq);
        alphapq=fangle(Op, Oq);

        %Toeplitz Structure of the extra-diagonal block
        C = [-Nq-(-Np):-1:-Nq-Np];
        R = [-Nq+(Np):1:Nq+(Np)];
        T = toeplitz(C,R);

        for iOperator=1:nTypeOfOperator
            this_Operator = ScalarTypeOfOperator(iOperator);
            this_weight = Weight(iOperator);
            switch this_Operator
                case -1, %Custom operator
                    Apq = Apq + this_weight*varargin{2}(Op, ap, Np, Oq, aq, Nq, k, varargin{3:end});
                %case 0 and 1 are the null matrix
                case 2, %Single Layer
                    Jm_kap = besselj(MNp,k*ap);
                    Jm_kaq = besselj(MNq,k*aq);
                    Apq = Apq + this_weight*1i*pi*0.5*sqrt(ap*aq)*diag(Jm_kap)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(Jm_kaq);
                case 3, % Double Layer
                    Jm_kap = besselj(MNp,k*ap);
                    dJm_kaq = dbesselj(MNq,k*aq);
                    Apq = Apq + this_weight*(-1i)*k*pi*0.5*sqrt(ap*aq)*diag(Jm_kap)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(dJm_kaq);
                case 4, % Dn Single Layer
                    dJm_kap = dbesselj(MNp,k*ap);
                    Jm_kaq = besselj(MNq,k*aq);
                    Apq = Apq + this_weight*1i*k*pi*0.5*sqrt(ap*aq)*diag(dJm_kap)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(Jm_kaq);
                case 5, % Dn Double Layer
                    dJm_kap = dbesselj(MNp,k*ap);
                    dJm_kaq = dbesselj(MNq,k*aq);
                    Apq = Apq + this_weight*(-1i)*k^2*pi*0.5*sqrt(ap*aq)*diag(dJm_kap)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(dJm_kaq);
                case 6, % Single scattering preconditioned integral operator for Dirichlet
                    H1m_kap_inv = 1./besselh(MNp,1,k*ap);
                    Jm_kaq = besselj(MNq,k*aq);
                    Apq = Apq + this_weight*sqrt(aq/ap)*diag(H1m_kap_inv)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(Jm_kaq);
                case 7, % Single scattering preconditioned integral operator for Neumann
                    dH1m_kap_inv = 1./dbesselh(MNp,1,k*ap);
                    dJm_kaq = dbesselj(MNq,k*aq);
                    Apq = Apq + this_weight*sqrt(aq/ap)*diag(dH1m_kap_inv)*besselh(T,1,k*bpq).*exp(1i.*T.*alphapq)*diag(dJm_kaq);
            end
        end
    end
end