function [CI_y,pred_inf,pred_sup,Hii] = ...
    LSQ_intervalo_prediccion(yPred,C,Z,sigma_y_hat,tcrit)

% Identifier: LSQ_intervalo_prediccion
%
% Input:
%   yPred       : valores predichos
%   C           : matriz de diseno o Jacobiano final
%   Z           : matriz de prediccion
%   sigma_y_hat : dispersion estandar estimada
%   tcrit       : valor critico de Student
%
% Output:
%   CI_y     : semiamplitud de prediccion
%   pred_inf : limite inferior
%   pred_sup : limite superior
%   Hii      : elementos diagonales de H

yPred = yPred(:);

[m,n] = size(C);

if m <= n || rank(C) < n
    error('La matriz C no cumple las condiciones.');
end

if size(Z,2) ~= n || size(Z,1) ~= numel(yPred)
    error('Dimensiones incompatibles.');
end

% Matriz de las ecuaciones normales
A = C.'*C;

% Inversa mediante resolucion de sistema lineal
M = A \ eye(n);

% Calculo de Hii
q = size(Z,1);

Hii = zeros(q,1);

for i = 1:q

    z = Z(i,:).';

    % Hii = z'*(C'*C)^(-1)*z
    Hii(i) = z.'*M*z;

end

% Semiamplitud del intervalo de prediccion
CI_y = sqrt(1 + Hii)*sigma_y_hat*tcrit;

% Limites del intervalo
pred_inf = yPred - CI_y;
pred_sup = yPred + CI_y;

end