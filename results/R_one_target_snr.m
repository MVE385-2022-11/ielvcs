%% Initialize
clear variables
simname = 'one_target_snr';

%% Setup all simulations to be run
% Row parameters
r_N = [130];
r_SNR = [40 20 10 5];
r_v = [7.5];
r_win = {@hamming @hann @welch};
    s_winname = {'Hamming', 'Hann', 'Welch'};
    
% With and without oversampling;
Nprime1 = 128;
Nprime2 = 1024;

% Lambda adaptable for each Np and window.
% Should be adjusted for SNR...
% Selected at 10 m/s (and double checked at 5 m/s)
lambdas = zeros(2,length(r_win),length(r_SNR));
                  %Ham %Han %Wel
lambdas(:,:,1) = [  3   5   5   % 
                    5   5   5]; % oversampling
lambdas(:,:,2) = [ 30  30  30   % 
                   30  30  30]; % oversampling
% Hamming is more or less a lost cause at this point for os.
lambdas(:,:,3) = [ 60  60  60   % 
                   60  60  60]; % oversampling
% Here the components drown in noise
lambdas(:,:,4) = [100  170 100   % 
                  100  170 100]; % oversampling

% Parameter grid
sim_grid = cell(1, 4);
[sim_grid{:}] = ndgrid(r_N, 1:length(r_win), 1:length(r_SNR), r_v);
sim_grid = reshape(cat(4,sim_grid{:}),[],4);

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
span = @(X) fill([min(X) min(X) max(X) max(X)], [min(ylim) max(ylim) max(ylim) min(ylim)], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

for row=1:R
    if ~mask(row)
        continue
    end
    
    timer = tic;
    
    gr = num2cell(sim_grid(row,:));
    [N, winix, snrix, v] = deal(gr{:});
    snr = r_SNR(snrix);

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
    s = smod.s(amp, v);
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
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(1,winix,snrix), G1, tol);
        Z_LF(:,i) = fft(BNpNp1*beta);
        [vests, ~, ~] = spectrum_peak_est(Z_LF(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix] = min(abs(vests - v));
            vest(i,4) = vests(minix);
        end
        
        % LASSO-FISTA (oversampled)
        beta0 = zeros(size(G2,2),1);
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(2,winix,snrix), G2, tol);
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
        s_winname{winix}, snr, v,...
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
        if ~any(isnan(vest(:,1)))
            xline(mean(vest(:,1),'omitnan'),'r');
            span(vest(:,1));
        end
        title(sprintf("SP: b = %.3f ± %.5f", b(row,1), b_err(row,1)))
    
        % MF
        nexttile
        dbplot_CI(Z_MF, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        Q = quantile(vest(:,2), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~any(isnan(vest(:,2)))
            xline(mean(vest(:,2),'omitnan'),'r');
            span(vest(:,2));
        end
        title(sprintf("MF: N' = %d, b = %.3f ± %.5f", Nprime1, b(row,2), b_err(row,2)))
        
        % MF (oversampled)
        nexttile
        dbplot_CI(Z_MF_O, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        Q = quantile(vest(:,3), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~any(isnan(vest(:,3)))
            xline(mean(vest(:,3),'omitnan'),'r');
            span(vest(:,3));
        end
        title(sprintf("MF: N' = %d, b = %.3f ± %.5f", Nprime2, b(row,3), b_err(row,3)))
    
        % Info
        nexttile
        axis ij
        set(gca,'Visible','off')
        text(0,0,sprintf("K = %d\nN = %d\nv = %.3f\nSNR = %.3f\nwindow = %s",...
            K, N, v, snr, s_winname{winix}));
        
        % LASSO-FISTA
        nexttile
        dbplot_ALL(Z_LF, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        Q = quantile(vest(:,4), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~any(isnan(vest(:,4)))
            xline(mean(vest(:,4),'omitnan'),'r');
            span(vest(:,4));
        end
        title(sprintf("LASSO: \x03bb = %g, b = %.3f ± %.5f", lambdas(1,winix,snrix), b(row,4), b_err(row,4)))
        
        % LASSO-FISTA (oversampled)
        nexttile
        dbplot_ALL(Z_LF_O, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        Q = quantile(vest(:,5), [0.025 0.5 0.975]);
        xline(v,'k');
        if ~any(isnan(vest(:,5)))
            xline(mean(vest(:,5),'omitnan'),'r');
            span(vest(:,5));
        end
        title(sprintf("LASSO: \x03bb = %g, b = %.3f ± %.5f", lambdas(2,winix,snrix), b(row,5), b_err(row,5)))

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
blims = [min(b_f - b_err,[],'all') max(b_f+b_err,[],'all')];

f = figure();
t = tiledlayout(f, 1, 3);
xlabel(t,'Target integrated SNR (dB)','FontSize',9);
ylabel(t,'Bias (m/s)','FontSize',9);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for winix = 1:3
    ax = nexttile;
    hold on
    yline(0, 'k', 'Alpha', 0.5);
    errorbar(r_SNR, b_f(winix:3:R,1), b_err(winix:3:R,1), '-o', 'Color', [217,95,2]./255);
    errorbar(r_SNR, b_f(winix:3:R,2), b_err(winix:3:R,2), '--o', 'Color', [27,158,119]./255);
    errorbar(r_SNR, b_f(winix:3:R,3), b_err(winix:3:R,3), '-o', 'Color', [27,158,119]./255);
    errorbar(r_SNR, b_f(winix:3:R,4), b_err(winix:3:R,4), '--o', 'Color', [117,112,179]./255);
    errorbar(r_SNR, b_f(winix:3:R,5), b_err(winix:3:R,5), '-o', 'Color', [117,112,179]./255);
    ylim(blims);
    xlim([min(r_SNR) max(r_SNR)]);
    title(s_winname{winix});
    grid on
end
l = legend(ax, {'', 'SP', 'MF', 'MF (os)', 'LASSO' 'LASSO (os)'}, 'Orientation','Horizontal');
l.Layout.Tile = 'North';

