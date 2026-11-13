function [d,OK,num_It,B,diagHist] = EVP_JacobiPC4(A,tol,maxIt)

[n,m] = size(A);

if n ~= m
    error('La matriz A debe ser cuadrada.');
end

% Verificar que A sea simetrica
symError = norm(A-A','fro');

if symError > 1e-12*max(1,norm(A,'fro'))
    error('El metodo de Jacobi requiere una matriz simetrica.');
end

% Matriz inicial
B = (A+A')/2;

OK = 0;
num_It = maxIt;

% ---------------------------------------------------------
% Historial de diagonales para grafica de convergencia
% ---------------------------------------------------------

diagHist = zeros(n,maxIt);

for k = 1:maxIt

    % -----------------------------------------------------
    % Norma de los elementos fuera de la diagonal
    % -----------------------------------------------------

    Boff = B - diag(diag(B));
    errorOff = norm(Boff,'fro');

    % Comprobar si ya converge
    if errorOff <= tol
        OK = 1;
        num_It = k-1;
        break
    end

    % -----------------------------------------------------
    % Buscar el mayor elemento fuera de la diagonal
    % -----------------------------------------------------

    M = abs(triu(B,1));

    [~,index] = max(M(:));

    [p,q] = ind2sub([n,n],index);

    % -----------------------------------------------------
    % Elementos del pivote
    % -----------------------------------------------------

    app = B(p,p);
    aqq = B(q,q);
    apq = B(p,q);

    % -----------------------------------------------------
    % Angulo de rotacion de Jacobi
    % -----------------------------------------------------

    theta = 0.5*atan2(2*apq,aqq-app);

    c = cos(theta);
    s = sin(theta);

    % -----------------------------------------------------
    % Matriz de rotacion de Jacobi
    % -----------------------------------------------------

    J = eye(n);

    J(p,p) = c;
    J(p,q) = s;
    J(q,p) = -s;
    J(q,q) = c;

    % -----------------------------------------------------
    % Transformacion de similaridad
    % -----------------------------------------------------

    B = J'*B*J;

    % Mantener simetria frente a errores de redondeo
    B = (B+B')/2;

    % El elemento seleccionado debe ser cero
    B(p,q) = 0;
    B(q,p) = 0;

    % -----------------------------------------------------
    % GUARDAR DIAGONAL DE ESTA ITERACION
    % -----------------------------------------------------

    diagHist(:,k) = diag(B);

    % -----------------------------------------------------
    % Comprobar convergencia despues de la rotacion
    % -----------------------------------------------------

    Boff = B - diag(diag(B));
    errorOff = norm(Boff,'fro');

    if errorOff <= tol
        OK = 1;
        num_It = k;
        break
    end

end

% ---------------------------------------------------------
% Recortar historial a las iteraciones realmente realizadas
% ---------------------------------------------------------

diagHist = diagHist(:,1:num_It);

% ---------------------------------------------------------
% Limpiar residuos numericos fuera de la diagonal
% ---------------------------------------------------------

mask = ~eye(n);

B(mask & abs(B) < tol) = 0;

% ---------------------------------------------------------
% Autovalores aproximados
% ---------------------------------------------------------

d = diag(B);

end