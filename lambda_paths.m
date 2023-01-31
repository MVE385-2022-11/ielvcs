function [lambdas,paths]=lambda_paths(y, smod, SP, Nprime, lambdaRange, tol)
SPy = SP.F * y;
B = fourier_matrix(smod.N,Nprime);
G = SP.F * B;
beta0 = zeros(Nprime,1);

lambdas = logspace(log10(lambdaRange(1)), log10(lambdaRange(2)), 100);
paths = zeros(Nprime,length(lambdas));

for i = 1:length(lambdas)
    [beta, ~] = LASSO_FISTA(SPy, beta0, lambdas(i), G, tol);
    paths(:,i) = abs(beta);
    fprintf('.');
end
fprintf('\n');

% fake lasso struct
FitInfo = struct();
FitInfo.MSE = 0; % checked to validate
FitInfo.PredictorNames = arrayfun(@(v) sprintf('v=%.1f', v),...
    fftshift(smod.vp*(0:(Nprime-1))/Nprime - smod.vp/2),...
    'UniformOutput', false);
FitInfo.Lambda = lambdas;
FitInfo.Alpha = 1; % LASSO
[ax,fig] = lassoPlot(paths, FitInfo, PlotType='Lambda', XScale='log');
ax.YScale = 'log';
fig.Children(2).YScale = 'log';
fig.Children(2).XLabel.String = sprintf("df (vres = %.4f, vres' = %.4f)", smod.vp/smod.N, smod.vp/Nprime);
end
