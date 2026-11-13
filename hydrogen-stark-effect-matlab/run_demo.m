% Verify and visualize the m = 0 hydrogen Stark model.
root = fileparts(mfilename('fullpath'));
sourceFolder = fullfile(root,'src');
addpath(sourceFolder);
cleanup = onCleanup(@() rmpath(sourceFolder));

nmax = 4;
hartreeToEv = 27.211386245988;
fields = linspace(0,4e8,61);
tolerance = 5e-9;
maxIterations = 10000;
[H0,Z,~,~,~] = Hydrogen_Hamiltonian(nmax,0);
referenceAtZero = sort(repelem(-1./(2*(1:nmax).^2),1:nmax)).';
assert(norm(sort(eig(H0)) - referenceAtZero,inf) < 1e-12, ...
    'Field-free energies differ from the analytic reference.');

energies = zeros(size(H0,1),numel(fields));
errors = zeros(size(fields));
for k = 1:numel(fields)
    H = H0 + fields(k)/5.14220675112e11*Z;
    [values,ok] = EVP_JacobiPC4(H,tolerance,maxIterations);
    assert(ok == 1,'Jacobi did not reach the requested tolerance.');
    values = sort(values);
    errors(k) = norm(values - sort(eig(H)),inf);
    energies(:,k) = hartreeToEv*values;
end
assert(max(errors) < 1e-8,'Jacobi and eig disagree.');
midpoint = (numel(fields)+1)/2;
fprintf('Largest Jacobi/eig difference: %.3e Hartree\n',max(errors));
fprintf('n = 2 splitting at 2e8 V/m: %.5f eV\n', ...
    energies(3,midpoint)-energies(2,midpoint));

fig = figure('Color','w');
hold on
groups = {2:3,4:6,7:10};
colors = lines(numel(groups));
for group = 1:numel(groups)
    indices = groups{group};
    plot(fields/1e8,energies(indices(1),:),'Color',colors(group,:), ...
        'LineWidth',1.6,'DisplayName',sprintf('n = %d',group+1));
    for index = indices(2:end)
        plot(fields/1e8,energies(index,:),'Color',colors(group,:), ...
            'LineWidth',1.6,'HandleVisibility','off');
    end
end
grid on
xlabel('Electric field (10^8 V/m)')
ylabel('Energy (eV)')
title('Hydrogen Stark splitting (m = 0, n_{max} = 4)')
legend('Location','best')
exportgraphics(fig,fullfile(root,'stark-splitting.png'),'Resolution',200);
clear cleanup
