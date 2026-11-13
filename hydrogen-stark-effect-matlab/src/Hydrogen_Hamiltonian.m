function [H,Z,D,basis,f] = Hydrogen_Hamiltonian(nmax,F)

% ---------------------------------------------------------
% Hydrogen_Hamiltonian
%
% Input:
%   nmax  : maximum principal quantum number
%   F     : electric field [V/m]
%
% Output:
%   H     : dimensionless Hamiltonian matrix
%   Z     : dimensionless dipole matrix <i|z/a0|j>
%   D     : unperturbed Hamiltonian matrix
%   basis : quantum numbers [n l m]
%   f     : dimensionless electric field F/F0
%
% The calculation uses the m = 0 Stark block.
% ---------------------------------------------------------

% Atomic unit of electric field
F0 = 5.14220675112e11;       % [V/m]

% Dimensionless electric field
f = F/F0;

% ---------------------------------------------------------
% Generate basis automatically
%
% B = {|n,l,m=0>}
%
% n = 1,...,nmax
% l = 0,...,n-1
% ---------------------------------------------------------

N = nmax*(nmax+1)/2;

basis = zeros(N,3);

k = 0;

for n = 1:nmax

    for l = 0:n-1

        k = k + 1;

        basis(k,:) = [n l 0];

    end

end

% ---------------------------------------------------------
% Unperturbed Hamiltonian
%
% E_n / Eh = -1/(2*n^2)
% ---------------------------------------------------------

D = zeros(N);

for i = 1:N

    n = basis(i,1);

    D(i,i) = -1/(2*n^2);

end

% ---------------------------------------------------------
% Dimensionless dipole matrix
%
% Zij = <ni li 0 | z/a0 | nj lj 0>
% ---------------------------------------------------------

Z = zeros(N);

for i = 1:N

    ni = basis(i,1);
    li = basis(i,2);

    for j = i+1:N

        nj = basis(j,1);
        lj = basis(j,2);

        % Stark selection rule:
        % Delta l = +-1
        %
        % Delta m = 0 is automatically satisfied
        % because all basis states have m = 0.

        if abs(li-lj) == 1

            % ---------------------------------------------
            % Radial integral
            %
            % Integral R_i(rho) R_j(rho) rho^3 drho
            % ---------------------------------------------

            Ir = integral( ...
                @(rho) Hydrogen_Radial(ni,li,rho) .* ...
                       Hydrogen_Radial(nj,lj,rho) .* ...
                       rho.^3, ...
                0,Inf, ...
                'RelTol',1e-11, ...
                'AbsTol',1e-13);

            % ---------------------------------------------
            % Angular integral
            %
            % <li,0 | cos(theta) | lj,0>
            % ---------------------------------------------

            Iang = Stark_Angular(li,lj);

            % ---------------------------------------------
            % Complete dipole matrix element
            % ---------------------------------------------

            Zij = Ir*Iang;

            Z(i,j) = Zij;

            % Real symmetric matrix
            Z(j,i) = Zij;

        end

    end

end

% ---------------------------------------------------------
% Final dimensionless Hamiltonian
%
% H = D + f Z
% ---------------------------------------------------------

H = D + f*Z;

end