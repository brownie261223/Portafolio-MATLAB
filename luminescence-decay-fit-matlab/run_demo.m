% Fit an exponential to six observations from the IF321 exercise.
root = fileparts(mfilename('fullpath'));
sourceFolder = fullfile(root,'src');
addpath(sourceFolder);
cleanup = onCleanup(@() rmpath(sourceFolder));

time = [0; 0.5; 1; 1.5; 2; 2.5];
intensity = [5; 3.20; 2.05; 1.30; 0.85; 0.52];
model = @(t,p) p(1).*exp(-p(2).*t);
jacobian = @(t,p) [exp(-p(2).*t), -t.*model(t,p)];

[p,iterations,ok] = LSQ_gauss_newton_actualizado( ...
    model,jacobian,time,intensity,[4.8; 0.7],2000);
assert(ok == 1,'Gauss-Newton did not reach the stopping criterion.');
fitted = model(time,p);
residuals = intensity-fitted;
sse = sum(residuals.^2);
assert(norm(p-[5.00163746; 0.89424931],inf) < 5e-5 && ...
    abs(sse-0.00049575637) < 1e-8,'Fit differs from the reference.');

C = jacobian(time,p);
[~,~,~,~,sigma,tCritical] = LSQ_incertidumbre_parametros( ...
    p,intensity,fitted,C,0.05);
timeGrid = linspace(min(time),max(time),300).';
predicted = model(timeGrid,p);
[~,lower,upper] = LSQ_intervalo_prediccion( ...
    predicted,C,jacobian(timeGrid,p),sigma,tCritical);

fprintf('Amplitude: %.6f a.u. | Decay rate: %.6f 1/ms\n',p(1),p(2));
fprintf('Iterations: %d | SSE: %.8f | RMSE: %.6f a.u.\n', ...
    iterations,sse,sqrt(mean(residuals.^2)));

fig = figure('Color','w');
tiledlayout(2,1);
nexttile
fill([timeGrid; flipud(timeGrid)],[lower; flipud(upper)], ...
    [0.74 0.86 0.97],'EdgeColor','none','DisplayName','95% prediction interval');
hold on
plot(timeGrid,predicted,'Color',[0.10 0.32 0.61],'LineWidth',1.8, ...
    'DisplayName','Exponential fit');
scatter(time,intensity,28,'k','filled','DisplayName','Coursework observations');
ylabel('Intensity (a.u.)')
title('Luminescence decay')
legend('Location','northeast')
grid on

nexttile
stem(time,residuals,'filled');
yline(0,':')
xlabel('Time (ms)')
ylabel('Residual (a.u.)')
grid on
exportgraphics(fig,fullfile(root,'luminescence-fit.png'),'Resolution',200);
clear cleanup
