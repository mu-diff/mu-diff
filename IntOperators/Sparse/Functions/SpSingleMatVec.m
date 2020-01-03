% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Sparse matrix vector product between a single matrix A and a vector X
% The matrix A is stored as a cell and has the special structure:
% A{1} = LeftPart  (off-diagonal blocks + diagonal blocks)
% A{2} = Toeplitz Part (off-diagonal blocks)
% A{3} = RightPart  (off-diagonal blocks)
% -------------------------------------------------------------
% The matrix vector product is done through a cross-correlation using xcorr
% Matlab function.
% -------------------------------------------------------------
%   Y = SpSingleMatVec(X, M_modes, A)
%
% OUTPUT ARGUMENT:
% ----------------
% Y [size(X)] : A*X
%
% INPUT ARGUMENTS (N_scat = number of obstacles):
% -----------------------------------------------
%
% X       [sum(2*M_modes+1), 1] : Vector
% M_Modes [N_scat, 1]           : Truncation index in the Fourier series of
%                                 obstacles
% A       cell(3,1)             : Integral operator (sparse storage)
%
%
% See also SpMatVec, SpBlockIntegralOperator, SpIntegralOperator, xcorr
%
function Y = SpSingleMatVec(X, M_modes, A)

N_scat = length(M_modes);
Y = zeros(size(X));
Sp = 0;
for p=1:N_scat
    Np = M_modes(p);
    Sq = 0;
    for q =1: N_scat
        Nq = M_modes(q);
        if(q==p) %Diagonal part of the matrix.
            Xp = X(Sp+1: Sp + 2*Np+1);
            Y(Sp+1: Sp + 2*Np+1) = Y(Sp+1: Sp + 2*Np+1) + A{1}(1:2*Np+1,p,p).*Xp ;
        else
            Xq = X(Sq+1: Sq+2*Nq+1);
            %Roughly speaking :
            %A{1}(:,p,q) * A{2}(:,p,q) * A{3}(:,p,q) *Xq
            %Where A{2} is Toeplitz, and A{1} and A{3} are diagonals.
            Zq = A{3}(1:2*Nq+1,p,q).*Xq;
            SizeRootVector = 2*Np + 2*Nq +1;
%            Zq_extended = [Zq;zeros(SizeMax-SizeRootVector,1)];
            cross_cor = xcorr(Zq, conj(A{2}(1:SizeRootVector,p,q)));
            Y(Sp+1: Sp + 2*Np+1) = Y(Sp+1: Sp + 2*Np+1) + ...
                A{1}(1:2*Np+1,p,q).*cross_cor(2*Nq+1:2*Nq+1+2*Np);
            clear Zq Zq_extended cross_cor;
        end
        Sq = Sq + 2*Nq+1;
    end
    Sp = Sp + 2*Np+1;
end
