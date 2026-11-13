function R = Hydrogen_Radial(n,l,rho)

% ---------------------------------------------------------
% Dimensionless hydrogen radial wave function
%
% rho = r/a0
%
% R_nl(rho)
% ---------------------------------------------------------

k = n-l-1;

alpha = 2*l+1;

x = 2*rho/n;

% Normalization constant
C = (2/n^2)*sqrt( ...
    factorial(n-l-1)/factorial(n+l) );

% Generalized Laguerre polynomial
L = Generalized_Laguerre(k,alpha,x);

% Dimensionless radial wave function
R = C .* exp(-rho/n) .* ...
    (2*rho/n).^l .* L;

end