function [p, numIt, OK] = LSQ_gauss_newton_actualizado(...
    fun_g, fun_Jg, xdata, ydata, gamma0, maxIt)

xdata = xdata(:);
ydata = ydata(:);
p = gamma0(:);

numIt = 0;
OK = 0;

Ndec = 6;
TolP = 0.5*10^(-Ndec);
TolD = 0.5*10^(-Ndec);

% Residuo modelo - datos, coherente con b = -J.'*e.
e = fun_g(xdata, p) - ydata;
D = e.'*e;

for k = 1:maxIt
    J = fun_Jg(xdata, p);
    A = J.'*J;
    b = -J.'*e;

    delta = A \ b;
    pNuevo = p + delta;

    eNuevo = fun_g(xdata, pNuevo) - ydata;
    DNuevo = eNuevo.'*eNuevo;

    errorP = norm(pNuevo - p, 2);
    errorD = abs(DNuevo - D);

    p = pNuevo;
    numIt = k;

    if errorP <= TolP
        OK = 1;
        return;
    end

    if errorD <= TolD
        OK = 1;
        return;
    end

    e = eNuevo;
    D = DNuevo;
end
end