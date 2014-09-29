% mu-diff - Copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Compute the Calderon projector matrix C given by
% C = [I/2+M, L ; D, I/2+N]
% The result is
% [gamma_D u] = [I/2+M, L]*[gamma_D u]
% [gamma_N u]   [D, I/2+N] [gamma_N u]
% where gamma_D u = Dirichlet trace of u on Gamma (boundary)
% and gamma_N u = normal.grad(u) on Gamma, where normal is directed
% outwardly to the obstacle
% -------------------------------------------------------------------------
%   A = CalderonProjector(O, a, M_modes, k)
%
% Output Arguments :
% -----------------
% A [4*sum(M_modes+1) x 4*sum(M_modes+1)] : Matrix of the system
%
% Input Arguments (N_scat = number of obstacles):
% -----------------
% O       [2 x N_scat]  : Coordinates of the center of the disk 1
% a       [1 x N_scat]  : Radius of disk 1
% M_Modes [1 x N_scat]  : Truncation index in the Fourier series of
%                         obstacles
% k       [1 x 1]       : Wavenumber in the vacuum
% 
% See also BlockCalderonProjector, IntegralOperator,
% SpBlockIntegralOperator, SpIntegralOperator

function A = CalderonProjector(O, a, M_modes, k, varargin)

    %Initialization
    N_scat = length(a);
    nvarargin = length(varargin);
    if(nvarargin >= 1)
       Weight = varargin{1};
    else
       Weight = ones(N_scat,N_scat);
    end
    if(~isscalar(Weight) && (size(Weight,1)~=N_scat || size(Weight,2)~=N_scat))
        error('Weight must be a N_scatxN_scat matrix');
    end
    Sum_modes = 2.*M_modes + 1;
    sizeA = 2*sum(Sum_modes);
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
        MNp_bis = 2*Np+1+[-Np:Np].';
        MNp_twice = [MNp; MNp_bis];
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
            MNq_bis = 2*Nq+1+[-Nq:Nq].';
            MNq_twice = [MNq; MNq_bis];
            %Compute the block matrix
            this_k = GetK(k,p,q);
            if(Weight(p,q) ~=0)
                A(Sp + MNp_twice +(Np+1), Sq + MNq_twice +(Nq+1)) = BlockCalderonProjector(Op, ap, Np, Oq, aq, Nq, this_k);
            end
            %Upgrade of the column-counter
            Sq = Sq + 2*(2*Nq+1);
        end
        %Upgrade of the row-counter
        Sp = Sp + 2*(2*Np+1); 
    end
end


%%

function this_k = GetK(k,p,q)
    if(~isscalar(k))
        if(p~=q)
            this_k = k(p);
        else
            this_k = k(p);
        end
    else
        this_k = k;
    end
end