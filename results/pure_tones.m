%% Initialize
clear
clf

%% Settings
N = 512+2;
f0 = 3e9; 
fp = 1e4;

% Number of simulations
K = 1;

smod = sigmod(N, f0, fp);
win{1} = @hamming;
win{2} = @hann;
win{3} = @welch;
win_names = ["Ham" "Han" "Wel"];
amp = 1;
vs = [2.5, 5, 10];

t = tiledlayout(length(vs), length(win));
%title(t, sprintf('Pure signal, varying velocity and window, v_{res} = %.2f', smod.vp/smod.N));
xlabel(t, 'Velocity (m/s)', 'Fontsize', 9);
ylabel(t, 'Amplitude (dB)', 'Fontsize', 9);
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% Run
for i = 1:length(vs)
    Y = smod.s(amp, vs(i)); 

    for j = 1:length(win)
        SP = SPsetup(smod,[-1 2 -1]/sqrt(6), win{j}(N-2)); 
        % Signal processing
        Z_SP = zeros(SP.M,K);
        % Here we are targeting a peak.
        % This may give the wrong result if the peak is not detectable.
        % It also underestimates bias for difficult situations since it assumes
        % a peak in fact exists.
        vest_SP = nan;
        
        y = Y(:,1);
        % SP
        z = SP.F*y;
        Z_SP(:,K) = z;
        [vests, ~, ~] = spectrum_peak_est(z,smod.vp,-20,3);
        if ~isempty(vests)
            % Get closest peak
            [~,minix] = min(abs(vests - vs(i)));
            vest_SP(K) = vests(minix);
        end

        nexttile;
        dbplot(Z_SP, smod.vp);
        xlim([-40 40])
        ylim([-80 40])
        xline(vs(i),'k');
        xline(mean(vest_SP),'r');
        xlabel('');
        ylabel('');
        
        b_SP = vest_SP - vs(i);
        title(sprintf('%s, v = %.1f, b = %.3f', win_names(j), vs(i), mean(b_SP)));

    end
end
