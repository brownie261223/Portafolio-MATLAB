function L = Generalized_Laguerre(k,alpha,x)

% ---------------------------------------------------------
% Generalized Laguerre polynomial
%
% L_k^alpha(x)
%
% Calculated using recurrence relations.
% ---------------------------------------------------------

if k == 0

    L = ones(size(x));

    return

end

if k == 1

    L = 1 + alpha - x;

    return

end

L0 = ones(size(x));

L1 = 1 + alpha - x;

for j = 2:k

    L = ( ...
        (2*j-1+alpha-x).*L1 ...
        - (j-1+alpha).*L0 ...
        )/j;

    L0 = L1;

    L1 = L;

end

end