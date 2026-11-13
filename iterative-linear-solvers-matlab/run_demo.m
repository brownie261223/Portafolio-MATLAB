% Compare the IF321 Jacobi and Gauss-Seidel implementations.
root = fileparts(mfilename('fullpath'));
sourceFolder = fullfile(root,'src');
addpath(sourceFolder);
cleanup = onCleanup(@() rmpath(sourceFolder));

A = [3 -1 0 -1 0 0 0 0;
    -1 3 -1 0 -1 0 0 0;
     0 -1 2 0 0 0 0 -1;
    -1 0 0 3 -1 -1 0 0;
     0 -1 0 -1 3 0 -1 0;
     0 0 0 -1 0 2 -1 0;
     0 0 0 0 -1 -1 3 -1;
     0 0 -1 0 0 0 -1 2];
b = [99.9; -.15; -.08; -.12; -.20; -.10; -.15; -.10];
x0 = 100*ones(8,1);
stepTolerance = 5e-5;
alphas = [1; 0.5; 0.2; 0.1; 0.05; 0.01];
results = zeros(numel(alphas),8);

for k = 1:numel(alphas)
    alpha = alphas(k);
    A(1,1) = alpha+2;
    b(1) = 100*alpha-0.10;
    [xJ,okJ,nJ,historyJ] = LSE_JacobiPC3(A,b,x0,stepTolerance,20000);
    [xGS,okGS,nGS,historyGS] = LSE_GaussSeidelPC3(A,b,x0,stepTolerance,20000);
    assert(okJ && okGS,'An iterative solver did not meet the step criterion.');
    exact = A\b;
    rhoJ = max(abs(eig(eye(8)-diag(diag(A))\A)));
    rhoGS = max(abs(eig(eye(8)-tril(A)\A)));
    results(k,:) = [nJ, nGS, rhoJ, rhoGS, norm(A*xJ-b), ...
        norm(A*xGS-b), norm(xJ-exact), norm(xGS-exact)];
    if k == 1
        baseHistoryJ = historyJ;
        baseHistoryGS = historyGS;
    end
end

assert(isequal(results([1 end],1:2),[247 135; 15982 8689]), ...
    'Iteration counts differ from the independent reference.');
summary = array2table([alphas results],'VariableNames', ...
    {'Alpha','JacobiSteps','GaussSeidelSteps','JacobiRadius','GaussSeidelRadius', ...
     'JacobiResidual','GaussSeidelResidual','JacobiError','GaussSeidelError'});
disp(summary)

fig = figure('Color','w');
tiledlayout(2,1);
nexttile
semilogy(baseHistoryJ,'LineWidth',1.6,'DisplayName','Jacobi');
hold on
semilogy(baseHistoryGS,'LineWidth',1.6,'DisplayName','Gauss-Seidel');
yline(stepTolerance,'--','DisplayName','Step tolerance');
xlabel('Iteration')
ylabel('Step 2-norm')
title('Convergence at alpha = 1')
grid on
legend('Location','best')

nexttile
loglog(alphas,results(:,1),'-o','LineWidth',1.6,'DisplayName','Jacobi');
hold on
loglog(alphas,results(:,2),'-s','LineWidth',1.6,'DisplayName','Gauss-Seidel');
xlabel('Alpha')
ylabel('Iterations')
title('Iteration count as alpha changes')
grid on
legend('Location','best')
exportgraphics(fig,fullfile(root,'solver-convergence.png'),'Resolution',200);
clear cleanup
