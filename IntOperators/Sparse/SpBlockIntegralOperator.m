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
% SPARSE VERSION
% -------------------------------------------------------------------------
% [LeftPart, MiddlePart, RightPart] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, k, TypeOfOperator, ...)
%
% OUTPUT ARGUMENTS:
% -----------------
% LeftPart   [2*Nmax+1, 1]       : Left part in the off-diag. blocks and
%                                  also the diagonal blocks of the matrix
% MiddlePart [2*(2*Nmax+1)-1, 1] : Middle part (Toeplitz part) of the 
%                                  off-diag. blocks
% RightPart  [2*Nmax+1, 1]       : Right part of the off-diag. blocks
%
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
% Nmax            [1 x 1]    : Reference number N (max of every other Np)
% k               [1 x 1]    : Wavenumber in the vacuum
% TypeOfOperator (See below) : Specifies the integral operator. 
%                              See Comon/Parser.m for correspondance. 
%
% TypeOfOperator acceptable size/type:
% ------------------------------------
% - SINGLE VALUE: SpBlockIntegralOperator(..., 2) or SpBlockIntegralOperator(..., {'L'})
%
% OPTION:
% -------
% [...] = SpBlockIntegralOperator(..., TypeOfOperator, Weight)
% Weight  [1 x 1]    : multiplicative constant to apply to the operator   
%
% EXAMPLE:
% --------
% SpBlockIntegralOperator(..., 2, 0.5) will produce 0.5*L
%   OR: SpBlockIntegralOperator(..., {'L'}, 0.5) produces 0.5*L
% 
% See also Parser, IntegralOperator,
% BlockIntegralOperator, SpIntegralOperator, BlockSingleLayer,
% BlockDnSingleLayer, BlockDoubleLayer, BlockDnDoubleLayer,
% BlockPrecondDirichlet, BlockPrecondNeumann

%%

function [LeftPart, MiddlePart, RightPart] = SpBlockIntegralOperator(Op, ap, Np, Oq, aq, Nq, Nmax, k, TypeOfOperator, varargin)

    ScalarTypeOfOperator = ParserIntegralOperator(TypeOfOperator);
    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(size(ScalarTypeOfOperator));
    end

    if(length(Weight) ~=length(ScalarTypeOfOperator))
        error('Matrix of TypeOfOperator and Weight must be of same size!');
    end
    if(length(ScalarTypeOfOperator) > 1)
        error('TypeOfOperator must be a scalar in SpBlockIntegralOperator!');
    end
    SizeMax = 2*Nmax+1;
    LeftPart = zeros(SizeMax,1);
    MiddlePart = zeros(2*SizeMax-1,1);
    RightPart = zeros(SizeMax,1);
    %Vector of the modes associated with the scatterer 1 and 2
    MNp = [-Np:Np].';

    if (Op(1) == Oq(1) && Op(2) == Oq(2) && ap == aq) % Diagonal block (which is diagonal)
        switch ScalarTypeOfOperator
          case -1,%Custom operator
                [LeftPart(1:2*Np+1), RightPart(1:2*Np+1)] = ...
                    varargin{2}(Op, ap, Np, Oq, aq, Nq, Nmax, k, varargin{3:end});            
            case 0, %Null matrix
            case 1, %Identity
                LeftPart(1:2*Np+1) = ones(2*Np+1,1);
            case 2, % Single Layer
                Jm_kap = besselj(MNp,k*ap);
                H1m_kap = besselh(MNp,1,k*ap);
                LeftPart(1:2*Np+1) = 1i*ap*pi*0.5*H1m_kap.*Jm_kap;
            case 3, % Double Layer
                dJm_kap = dbesselj(MNp,k*ap);
                H1m_kap = besselh(MNp,1,k*ap);
                LeftPart(1:2*Np+1) = 0.5 - 1i*k*ap*pi*0.5*H1m_kap.*dJm_kap;
            case 4, % Dn Single Layer
                dJm_kap = dbesselj(MNp,k*ap);
                H1m_kap = besselh(MNp,1,k*ap);
                LeftPart(1:2*Np+1) = -0.5 + 1i*k*ap*pi*0.5*H1m_kap.*dJm_kap;
            case 5, % Dn Double Layer
                dJm_kap = dbesselj(MNp,k*ap);
                dH1m_kap = dbesselh(MNp,1,k*ap);
                LeftPart(1:2*Np+1) = -1i*k^2*ap*pi*0.5*dH1m_kap.*dJm_kap;
            case 6, % Single scattering preconditioned single layer (Dirichlet)
                LeftPart(1:2*Np+1) = ones(2*Np+1,1);
            case 7, % Single scattering preconditioned dn double layer (Neumann)
                LeftPart(1:2*Np+1) = ones(2*Np+1,1);
        end
        %Weight is carried out by the left part
        if(Weight ~= 1)
            LeftPart(1:2*Np+1) = Weight*LeftPart(1:2*Np+1);
        end
    else %Off-diagonal blocks
        MNq = [-Nq:Nq].';
        % Common part to every integral operators:
        % Distance and angle between centers
        bpq = norm(Op - Oq);
        alphapq=fangle(Op, Oq);
        %Toeplitz Structure of the extra-diagonal block
        C = [-Nq+Np:-1:-Nq-Np];
        R = [-Nq+Np:1:Nq+Np];
    %    [Nq-(-Np):-1:-Nq-(-Np), -Nq+(Np)-1:-1:-Nq-(Np)];
    %Like a convolution product
    %    Root_Vector = [R(end:-1:1), C(2:end)].';
    %or a cross-corellation
        Root_Vector = [C(end:-1:2), R].';
        if(ScalarTypeOfOperator > 1)
            MiddlePart(1:length(Root_Vector)) = 1i*pi*0.5*besselh(Root_Vector,1,k*bpq).*exp(1i.*Root_Vector.*alphapq);
        else
            MiddlePart(1:length(Root_Vector)) = zeros(size(Root_Vector));
        end
        %Uncommon part:
        switch ScalarTypeOfOperator
          case -1, %custom operator
            [LeftPart(1:2*Np+1), RightPart(1:2*Nq+1)] = ...
                    varargin{2}(Op, ap, Np, Oq, aq, Nq, Nmax, k, varargin{3:end});
            case 0, %Null matrix
            case 1, %Identity (=> null on off diagonal blocks)
            case 2, %Single Layer
                Jm_kap = besselj(MNp,k*ap);
                Jm_kaq = besselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = sqrt(ap)*Jm_kap;
                RightPart(1:2*Nq+1) = sqrt(aq)*Jm_kaq;
            case 3, %Double Layer
                Jm_kap = besselj(MNp,k*ap);
                dJm_kaq = dbesselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = -sqrt(ap)*Jm_kap;
                RightPart(1:2*Nq+1) = k*sqrt(aq)*dJm_kaq;
            case 4, %Dn Single Layer
                dJm_kap = dbesselj(MNp,k*ap);
                Jm_kaq = besselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = k*sqrt(ap)*dJm_kap;
                RightPart(1:2*Nq+1) = sqrt(aq)*Jm_kaq;
            case 5, %Dn Double Layer
                dJm_kap = dbesselj(MNp,k*ap);
                dJm_kaq = dbesselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = -k*sqrt(ap)*dJm_kap;
                RightPart(1:2*Nq+1) = k*sqrt(aq)*dJm_kaq;
            case 6, % Single scattering preconditioned single layer (Dirichlet)
                H1m_kap_inv = 1./besselh(MNp,1,k*ap);
                Jm_kaq = besselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = 1./sqrt(ap)*H1m_kap_inv;
                RightPart(1:2*Nq+1) = 1/(1i*pi*0.5)*sqrt(aq)*Jm_kaq;
            case 7, % Single scattering preconditioned dn double layer (Neumann)
                dH1m_kap_inv = 1./dbesselh(MNp,1,k*ap);
                dJm_kaq = dbesselj(MNq,k*aq);
                LeftPart(1:2*Np+1) = 1./sqrt(ap)*dH1m_kap_inv;
                RightPart(1:2*Nq+1) = 1/(1i*pi*0.5)*sqrt(aq)*dJm_kaq;
        end
        %Weight is carried out by the left part
        if(Weight ~=1)
            LeftPart(1:2*Np+1) = Weight*LeftPart(1:2*Np+1);
        end
    end
end
