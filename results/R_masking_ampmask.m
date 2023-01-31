%% Initialize
clear variables
simname = 'masking_ampmask_20';

%% Setup all simulations to be run
% Row parameters
r_N = [130];
r_SNR = [20];
r_v = [10];
r_win = {@hamming @hann @welch};
    s_winname = {'Hamming', 'Hann', 'Welch'};
r_ampmask = 10.^[0 1 1.5 2];

v_mask = 100;
    
% With and without oversampling;
Nprime1 = 128;
Nprime2 = 1024;

% Consider recalibrating lambda
% Also should perhaps adapt to ampmask (but is that realistic?)
lambdas = zeros(2,length(r_win),length(r_ampmask));
                  %Ham %Han %Wel
%lambdas(:,:,1) = [ 3    5    5    % 
%                   5    5    5];  % oversampling
%lambdas(:,:,2) = [ 3    5    5    % 
%                   5    5    5];  % oversampling
%lambdas(:,:,3) = [ 10    10   10  % 
%                   20    5    5]; % oversampling
%lambdas(:,:,4) = [ 30   30  30    % 
%                   20   20  10];  % oversampling
lambdas(:,:,1) = [ 30 30 30    % 
                   30 30 30];  % oversampling
lambdas(:,:,2) = [ 30 30 30    % 
                   30 30 30];  % oversampling
lambdas(:,:,3) = [ 30 30 30  % 
                   30 30 30]; % oversampling
lambdas(:,:,4) = [ 30 30 30    % 
                   30 30 30];  % oversampling

% Parameter grid
sim_grid = cell(1, 5);
[sim_grid{:}] = ndgrid(r_N, 1:length(r_win), 1:length(r_SNR), r_v, 1:length(r_ampmask));
sim_grid = reshape(cat(5,sim_grid{:}),[],5);

%% Settings
% General params
f0 = 3e9; 
fp = 1e4;

K = 150; % Number of simulations
tol = 1e-4; % FISTA tolerance
ip = 3; % Interpolation pts

% Output file
if ~exist('out','dir')
    mkdir('out');
    mkdir('out/fig');
end
f = fopen(['out/' simname '.log'], 'a');
fprintf(f, '%%\n%% %s Simulation %s\n', simname, datetime);
fprintf(f, "%% K = %d; N = [%s]; N' = [%d %d]; lambdas = [%s]\n",...
        K, sprintf('%d ', r_N), Nprime1, Nprime2, sprintf('%g ', lambdas));

% Options
R = size(sim_grid,1);
% Set entry to zero to skip running.
mask = ones(R, 1);
% Seed. For reproducibility. Here we use same noise pattern for every
% setting to avoid random effects across speeds/windows.
%- rng(123);
%- seed = randi(1000,R,1);
seed = 123*ones(R,1);
% MinPeakHeight in dB. To avoid spurious identifications
mph = -20;

%% Run
b = zeros(R,5);
b_err = zeros(R,5);
b_sucfrac = zeros(R,5);

% function used in some plots
span = @(X) fill([min(X,[],'omitnan') min(X,[],'omitnan') max(X,[],'omitnan') max(X,[],'omitnan')], [min(ylim) max(ylim) max(ylim) min(ylim)], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

for row=1:R
    if ~mask(row)
        continue
    end
    
    timer = tic;
    
    gr = num2cell(sim_grid(row,:));
    [N, winix, snrix, v, ampmaskix] = deal(gr{:});
    snr = r_SNR(snrix);
    ampmask = r_ampmask(ampmaskix);

    smod = sigmod(N, f0, fp);
    SP = SPsetup(smod, [-1 2 -1]/sqrt(6), r_win{winix}(N-2));

    amp = 1;
    sigma = noise_isnr(snr, N, amp);
    
    % Processing preparation
    BNNp1 = fourier_matrix(smod.N,Nprime1);
    BNpNp1 = dftmtx(Nprime1)';
    G1 = SP.F * BNNp1;
    BNNp2 = fourier_matrix(smod.N,Nprime2);
    BNpNp2 = dftmtx(Nprime2)';
    G2 = SP.F * BNNp2;
  
    Z_SP   = zeros(SP.M, K);
    Z_MF   = zeros(size(G1,2), K);
    Z_LF   = zeros(size(G1,2), K);
    Z_MF_O = zeros(size(G2,2), K);
    Z_LF_O = zeros(size(G2,2), K);

    % Simulate K samples
    rng(seed(row));
    if v ~= 0
        s = smod.s([amp ampmask], [v v_mask]);
    else
        s = smod.s(ampmask, v_mask);
    end
    if sigma == 0
        n = zeros(smod.N,K);
    else
        n = smod.randcn(sigma, smod.N, K);
    end
    Y = s + n;
    
    % Prepare peak estimation array
    % We will use the fact that we know the true location.
    % Even the smallest peak will be identified. We do not try to
    % simulate distinguishability over noise.
    vest = ones(K,5)*nan;

    for i=1:K
        y = Y(:,i);
    
        % SP
        z = SP.F*y;
        Z_SP(:,i) = z;
        [vests, ~, ~] = spectrum_peak_est(z,smod.vp,mph,ip);
        if ~isempty(vests)
            % Get closest peak
            [~,minix] = min(abs(vests - v));
            vest(i,1) = vests(minix);
        end
    
        % MF
        Z_MF(:,i) = G1'*z;
        [vests, ~, ~] = spectrum_peak_est(Z_MF(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix] = min(abs(vests - v));
            vest(i,2) = vests(minix);
        end
        
        % MF oversampled
        Z_MF_O(:,i) = G2'*z;
        [vests, ~, ~] = spectrum_peak_est(Z_MF_O(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix] = min(abs(vests - v));
            vest(i,3) = vests(minix);
        end
    
        % LASSO-FISTA
        beta0 = zeros(size(G1,2),1);
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(1,winix,ampmaskix), G1, tol);
        Z_LF(:,i) = fft(BNpNp1*beta);
        [vests, ~, ~] = spectrum_peak_est(Z_LF(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix] = min(abs(vests - v));
            vest(i,4) = vests(minix);
        end
        
        % LASSO-FISTA (oversampled)
        beta0 = zeros(size(G2,2),1);
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(2,winix,ampmaskix), G2, tol);
        Z_LF_O(:,i) = fft(BNpNp2*beta);
        [vests, ~, ~] = spectrum_peak_est(Z_LF_O(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix] = min(abs(vests - v));
            vest(i,5) = vests(minix);
        end
    end

    b(row,:) = mean(vest - v, 1, 'omitnan');
    b_sucfrac(row,:) = 1 - mean(isnan(vest), 1);
    b_err(row,:) = std(vest - v, 1, 'omitnan')./sqrt(b_sucfrac(row,:)*K);

    fprintf(f, "%s & %f & %f & \n\t %f ± %g (%.2f) & \n\t %f ± %g (%.2f) & %f ± %g (%.2f) & \n\t %f ± %g (%.2f) & %f ± %g (%.2f) \\\\\n",...
        s_winname{winix}, ampmask, v,...
        reshape([b(row,:); b_err(row,:); b_sucfrac(row,:)], 1, 15));
    
    fprintf('%d/%d in %f\n', row, R, toc(timer));
    
        for fig_section = 1
        % Figure
        fig = figure('Position', [0 0 1440 900], 'WindowState', 'minimized');
        set(fig,'defaulttextInterpreter','none')
        t = tiledlayout(fig,2,3);
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
        xlabel(t, 'Velocity (m/s)','Fontsize',9);
        ylabel(t, 'Amplitude (dB)','Fontsize',9);

        % SP
        nexttile
        dbplot_CI(Z_SP, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
    
        Q = quantile(vest(:,1), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~all(isnan(vest(:,1)))
            xline(mean(vest(:,1),'omitnan'),'r');
            span(vest(:,1));
        end
        title(sprintf("SP: b = %.3f ± %.5f (%.2f)", b(row,1), b_err(row,1), b_sucfrac(row,1)));
    
        % MF
        nexttile
        dbplot_CI(Z_MF, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        Q = quantile(vest(:,2), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~all(isnan(vest(:,2)))
            xline(mean(vest(:,2),'omitnan'),'r');
            span(vest(:,2));
        end
        title(sprintf("MF: b = %.3f ± %.5f (%.2f)", b(row,2), b_err(row,2), b_sucfrac(row,2)));
        
        % MF (oversampled)
        nexttile
        dbplot_CI(Z_MF_O, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        Q = quantile(vest(:,3), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~all(isnan(vest(:,3)))
            xline(mean(vest(:,3),'omitnan'),'r');
            span(vest(:,3));
        end
        title(sprintf("MF (os): b = %.3f ± %.5f (%.2f)", b(row,3), b_err(row,3), b_sucfrac(row,3)));
    
        % Info
        nexttile
        axis ij
        set(gca,'Visible','off')
        text(0,0.5,sprintf("K = %d\nN = %d\nN' = [%s]\nv = %.3f\nSNR = %.3f\nwindow = %s\nampmask = %.3f\nvmask = %.3f",...
            K, N, sprintf("%d, %d", Nprime1, Nprime2), v, snr, s_winname{winix}, ampmask, v_mask), 'FontSize', 9);
        
        % LASSO-FISTA
        nexttile
        dbplot_ALL(Z_LF, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        Q = quantile(vest(:,4), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~all(isnan(vest(:,4)))
            xline(mean(vest(:,4),'omitnan'),'r');
            span(vest(:,4));
        end
        title(sprintf("LASSO: \x03bb = %g, b = %.3f ± %.5f (%.2f)", lambdas(1,winix), b(row,4), b_err(row,4), b_sucfrac(row,4)));
        
        % LASSO-FISTA (oversampled)
        nexttile
        dbplot_ALL(Z_LF_O, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        Q = quantile(vest(:,5), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~all(isnan(vest(:,5)))
            xline(mean(vest(:,5),'omitnan'),'r');
            span(vest(:,5));
        end
        title(sprintf("LASSO (os): \x03bb = %g, b = %.3f ± %.5f (%.2f)", lambdas(2,winix), b(row,5), b_err(row,5), b_sucfrac(row,5)));

        savefig(fig,sprintf('out/fig/%s.%d.fig', simname, row),'compact');
        close(fig);
    end
end

fclose(f);
save(['out/' simname '.mat']);

% Celebrate success
endsound_spec = {'y','Fs'};
endsound = load('train.mat',endsound_spec{:});
sound(endsound.y, endsound.Fs);

%% Plot lines
close all
nanfilt = [0 nan];
b_f = b + nanfilt((b_sucfrac < 0.95) + 1);
b_f(abs(b_f) > 15)= nan;
blims = [min(b_f - b_err,[],'all') max(b_f+b_err,[],'all')];

f = figure();
t = tiledlayout(f, 1, 3);
xlabel(t,'Masking object relative amplitude (dB)','FontSize',9);
ylabel(t,'Bias (m/s)','FontSize',9);
t.TileSpacing = 'compact';
t.Padding = 'compact';
r_ampmask_dB = 20*log10(r_ampmask);
for winix = 1:3
    ax = nexttile;
    hold on
    yline(0, 'k', 'Alpha', 0.5);
    errorbar(r_ampmask_dB, b_f(winix:3:R,1), b_err(winix:3:R,1), '-o', 'Color', [217,95,2]./255);
    errorbar(r_ampmask_dB, b_f(winix:3:R,2), b_err(winix:3:R,2), '--o', 'Color', [27,158,119]./255);
    errorbar(r_ampmask_dB, b_f(winix:3:R,3), b_err(winix:3:R,3), '-o', 'Color', [27,158,119]./255);
    errorbar(r_ampmask_dB, b_f(winix:3:R,4), b_err(winix:3:R,4), '--o', 'Color', [117,112,179]./255);
    errorbar(r_ampmask_dB, b_f(winix:3:R,5), b_err(winix:3:R,5), '-o', 'Color', [117,112,179]./255);
    ylim(blims);
    title(s_winname{winix});
    grid on
end
l = legend(ax, {'', 'SP', 'MF', 'MF (os)', 'LASSO' 'LASSO (os)'}, 'Orientation','Horizontal');
l.Layout.Tile = 'North';
