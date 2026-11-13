function [x, OK, NumIt, errorHist] = LSE_JacobiPC3(A, b, x0, Tol_a, MaxIt)

n = size(A,1);

OK = 0;
NumIt = MaxIt;

x = zeros(n,1);

errorHist = zeros(MaxIt,1);

for k = 1:MaxIt

    for i = 1:n

        suma1 = 0;
        suma2 = 0;

        for j = 1:i-1
            suma1 = suma1 + A(i,j)*x0(j);
        end

        for j = i+1:n
            suma2 = suma2 + A(i,j)*x0(j);
        end

        x(i) = (b(i) - suma1 - suma2)/A(i,i);

    end

    suma = 0;

    for i = 1:n
        suma = suma + (x(i) - x0(i))^2;
    end

    error = sqrt(suma);

    errorHist(k) = error;

    if error <= Tol_a
        OK = 1;
        NumIt = k;
        break
    end
    x0 = x;
end

errorHist = errorHist(1:NumIt);

end