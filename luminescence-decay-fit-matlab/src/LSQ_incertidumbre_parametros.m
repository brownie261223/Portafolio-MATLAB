function [u_p,CI_p,p_inf,p_sup,sigma_y_hat,tcrit] = ...
    LSQ_incertidumbre_parametros(p,ydata,yfit,C,alphaConf)

% Identifier: LSQ_incertidumbre_parametros
%
% Input:
%   p         : parametros ajustados
%   ydata     : datos experimentales
%   yfit      : valores ajustados
%   C         : matriz de diseno o Jacobiano final
%   alphaConf : nivel de significancia
%
% Output:
%   u_p         : incertidumbres estandar
%   CI_p        : semiamplitudes de confianza
%   p_inf       : limites inferiores
%   p_sup       : limites superiores
%   sigma_y_hat : dispersion estandar estimada
%   tcrit       : valor critico de Student

p = p(:);
ydata = ydata(:);
yfit = yfit(:);

[m,n] = size(C);
nu = m - n;

if nu <= 0 || rank(C) < n
    error('Se requiere m > n y rango completo.');
end

if numel(p) ~= n || numel(ydata) ~= m || numel(yfit) ~= m
    error('Dimensiones incompatibles.');
end

if alphaConf <= 0 || alphaConf >= 1
    error('alphaConf debe pertenecer a (0,1).');
end

% Residuos y dispersion
e = ydata - yfit;
D = e.'*e;

sigma_y_hat = sqrt(D/nu);

% Valor critico de Student
beta = betaincinv(alphaConf,nu/2,0.5);
tcrit = sqrt(nu*(1-beta)/beta);

% Matriz de covarianza de los parametros
A = C.'*C;
M = A \ eye(n);

cov_p = sigma_y_hat^2*M;

% Incertidumbres estandar
u_p = sqrt(diag(cov_p));

% Semiamplitudes de los intervalos
CI_p = tcrit*u_p;

% Limites de confianza
p_inf = p - CI_p;
p_sup = p + CI_p;

end