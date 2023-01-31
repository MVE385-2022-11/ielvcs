%% Initialize
clear variables
simname = 'separation_zoom';

%% Setup all simulations to be run
% Row parameters
r_N = [130];
r_SNR = [20];
r_v1 = linspace(5, 2.5, 4);
r_dv = linspace(5, 2.5, 4);
r_win = {@hann @welch};
    s_winname = {'Hann', 'Welch'};
    
% With and without oversampling;
Nprime1 = 128;
Nprime2 = 1024;

% Consider recalibrating lambda
% Todo reconfirm
lambdas = zeros(2,length(r_win),length(r_SNR));
                  %Han %Wel
lambdas(:,:,1) = [ 30   30   % 
                   30   30]; % oversampling

% Parameter grid
sim_grid = cell(1, 5);
[sim_grid{:}] = ndgrid(r_N, 1:length(r_win), 1:length(r_SNR), r_v1, r_dv);
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
sucfrac = zeros(R,5);
sucfrac_err = zeros(R,5);
b = zeros(R,5);


% function used in some plots
span = @(X) fill([min(X) min(X) max(X) max(X)], [min(ylim) max(ylim) max(ylim) min(ylim)], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

for row=1:R
    if ~mask(row)
        continue
    end
    
    timer = tic;
    
    gr = num2cell(sim_grid(row,:));
    [N, winix, snrix, v1, dv] = deal(gr{:});
    snr = r_SNR(snrix);
    v2 = v1 + dv;

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
    s = smod.s([amp amp], [v1 v2]);
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
    vest = ones(K,5,2)*nan;

    for i=1:K
        y = Y(:,i);
    
        % SP
        z = SP.F*y;
        Z_SP(:,i) = z;
        [vests, ~, ~] = spectrum_peak_est(z,smod.vp,mph,ip);
        if ~isempty(vests)
            % Get closest peaks
            [~,minix1] = min(abs(vests - v1));
            [~,minix2] = min(abs(vests - v2));
            vest(i,1,1) = vests(minix1);
            vest(i,1,2) = vests(minix2);
        end
    
        % MF
        Z_MF(:,i) = G1'*z;
        [vests, ~, ~] = spectrum_peak_est(Z_MF(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix1] = min(abs(vests - v1));
            [~,minix2] = min(abs(vests - v2));
            vest(i,2,1) = vests(minix1);
            vest(i,2,2) = vests(minix2);
        end
        
        % MF oversampled
        Z_MF_O(:,i) = G2'*z;
        [vests, ~, ~] = spectrum_peak_est(Z_MF_O(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix1] = min(abs(vests - v1));
            [~,minix2] = min(abs(vests - v2));
            vest(i,3,1) = vests(minix1);
            vest(i,3,2) = vests(minix2);
        end
    
        % LASSO-FISTA
        beta0 = zeros(size(G1,2),1);
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(1,winix,snrix), G1, tol);
        Z_LF(:,i) = fft(BNpNp1*beta);
        [vests, ~, ~] = spectrum_peak_est(Z_LF(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix1] = min(abs(vests - v1));
            [~,minix2] = min(abs(vests - v2));
            vest(i,4,1) = vests(minix1);
            vest(i,4,2) = vests(minix2);
        end
        
        % LASSO-FISTA (oversampled)
        beta0 = zeros(size(G2,2),1);
        [beta, ~] = LASSO_FISTA(z, beta0, lambdas(2,winix,snrix), G2, tol);
        Z_LF_O(:,i) = fft(BNpNp2*beta);
        [vests, ~, ~] = spectrum_peak_est(Z_LF_O(:,i),smod.vp,mph,ip);
        if ~isempty(vests)
            [~,minix1] = min(abs(vests - v1));
            [~,minix2] = min(abs(vests - v2));
            vest(i,5,1) = vests(minix1);
            vest(i,5,2) = vests(minix2);
        end
    end

    % Remove cases where the two peaks have completely merged.
    vest(abs(vest(:,:,1) - vest(:,:,2)).*ones(size(vest)) < 1e-5) = nan;
    sucfrac(row,:) = 1 - mean(isnan(vest(:,:,1)), 1);
    sucfrac_err(row,:) = sqrt(sucfrac(row,:) .* (1 - sucfrac(row,:))./K);
    % We still log estimated separation 
    b(row,:) = mean(vest(:,:,2) - vest(:,:,1) - (v2 - v1), 1, 'omitnan');

    fprintf(f, "%s & %f & %f & %f & \n\t %f ± %g (%.2f) & \n\t %.2f ± %g (%f) & %.2f ± %g (%f) & \n\t %.2f ± %g (%f) & %.2f ± %g (%f) \\\\\n",...
        s_winname{winix}, snr, v1, v2,...
        reshape([sucfrac(row,:); sucfrac_err(row,:); b(row,:)], 1, 15));
    
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
    
        xline(v1,'k');
        xline(v2,'k');
        if ~any(isnan(vest(:,1,:)))
            xline(mean(vest(:,1,1),'omitnan'),'r');
            xline(mean(vest(:,1,2),'omitnan'),'r');
            span(vest(:,1,1));
            span(vest(:,1,2));
        end
        title(sprintf("SP: %% = %.3f ± %.5f", sucfrac(row,1), sucfrac_err(row,1)))
    
        % MF
        nexttile
        dbplot_CI(Z_MF, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        xline(v1,'k');
        xline(v2,'k');
        if ~any(isnan(vest(:,2,:)))
            xline(mean(vest(:,2,1),'omitnan'),'r');
            xline(mean(vest(:,2,2),'omitnan'),'r');
            span(vest(:,2,1));
            span(vest(:,2,2));
        end
        title(sprintf("MF: N' = %d, %% = %.3f ± %.5f", Nprime1, sucfrac(row,2), sucfrac_err(row,2)))
        
        % MF (oversampled)
        nexttile
        dbplot_CI(Z_MF_O, smod.vp);
        xlim([-25 125]);
        ylim([-40 80]);
        xlabel('');
        ylabel('');
        
        xline(v1,'k');
        xline(v2,'k');
        if ~any(isnan(vest(:,3,:)))
            xline(mean(vest(:,3,1),'omitnan'),'r');
            xline(mean(vest(:,3,2),'omitnan'),'r');
            span(vest(:,3,1));
            span(vest(:,3,2));
        end
        title(sprintf("MF: N' = %d, %% = %.3f ± %.5f", Nprime2, sucfrac(row,3), sucfrac_err(row,3)))
    
        % Info
        nexttile
        axis ij
        set(gca,'Visible','off')
        text(0,0,sprintf("K = %d\nN = %d\nv1 = %.3f\nv2 = %.3f\ndv = %.3f\nSNR = %.3f\nwindow = %s",...
            K, N, v1, v2, dv, snr, s_winname{winix}));
        
        % LASSO-FISTA
        nexttile
        dbplot_ALL(Z_LF, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        xline(v1,'k');
        xline(v2,'k');
        if ~any(isnan(vest(:,4,:)))
            xline(mean(vest(:,4,1),'omitnan'),'r');
            xline(mean(vest(:,4,2),'omitnan'),'r');
            span(vest(:,4,1));
            span(vest(:,4,2));
        end
        title(sprintf("LASSO: \x03bb = %g, %% = %.3f ± %.5f", lambdas(1,winix), sucfrac(row,4), sucfrac_err(row,4)))
        
        % LASSO-FISTA (oversampled)
        nexttile
        dbplot_ALL(Z_LF_O, smod.vp);
        xlim([-25 125]);
        ylim([-100 max(ylim())]);
        xlabel('');
        ylabel('');
    
        xline(v1,'k');
        xline(v2,'k');
        if ~any(isnan(vest(:,5,:)))
            xline(mean(vest(:,5,1),'omitnan'),'r');
            xline(mean(vest(:,5,2),'omitnan'),'r');
            span(vest(:,5,1));
            span(vest(:,5,2));
        end
        title(sprintf("LASSO: \x03bb = %g, %% = %.3f ± %.5f", lambdas(2,winix), sucfrac(row,5), sucfrac_err(row,5)))

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
snrix = 1;
algname = {'SP', 'MF', 'MF (os)', 'LASSO', 'LASSO (os)'};

close all

f = figure();
t = tiledlayout(f, 3, 2);
xlabel(t,'Target 1 base velocity (m/s)','FontSize',9);
ylabel(t,'Target velocity difference (m/s)','FontSize',9);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for algix = [1 3 5]
for winix = 1:2
    baseix = winix;
    XX = baseix:2:R;
    ax = nexttile;
    C = reshape(sucfrac(XX,algix), length(r_v1), length(r_dv)).';
    
    %imagesc(r_v1, r_dv, C, [0 1]);
    pcolor(r_v1,r_dv,C);
    caxis([0 1]);
    shading('interp');
    hold on
    scatter(sim_grid(XX,4),sim_grid(XX,5),[],sucfrac(XX,algix),'filled','s',...
        'MarkerEdgeColor',[0 0 0], 'SizeData', 100);
    
    title(sprintf('%s %s', algname{algix}, s_winname{winix}));
end
end
cbh = colorbar; 
cbh.Layout.Tile = 'east';
