%% Initialize
clear variables
simname = 'interpolation_points';

%% Setup all simulations to be run
% Row parameters
r_N = [130 514];
r_v = [10 5 2.5];
r_win = {@hamming @hann @welch};
    s_winname = {'Hamming', 'Hann', 'Welch'};

c_ip = [3 5 7];

% Parameter grid
sim_grid = cell(1, 3);
[sim_grid{:}] = ndgrid(r_N, 1:length(r_win), r_v);
sim_grid = reshape(cat(3,sim_grid{:}),[],3);

%% Settings
% General params
f0 = 3e9; 
fp = 1e4;

% Output file
if ~exist('out','dir')
    mkdir('out');
    mkdir('out/fig');
end
f = fopen(['out/' simname '.log'], 'a');
fprintf(f, '%%\n%% %s Simulation %s\n', simname, datetime);

% Options
R = size(sim_grid,1);
% Set entry to zero to skip running.
mask = ones(R, 1);
% MinPeakHeight in dB. To avoid spurious identifications
mph = -20;

%% Run
b = zeros(R,length(c_ip));
b_err = zeros(R,length(c_ip));
b_sucfrac = zeros(R,length(c_ip));

% function used in some plots
span = @(X) fill([min(X) min(X) max(X) max(X)], [min(ylim) max(ylim) max(ylim) min(ylim)], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

for row=1:R
    if ~mask(row)
        continue
    end
    
    timer = tic;
    
    gr = num2cell(sim_grid(row,:));
    [N, winix, v] = deal(gr{:});

    smod = sigmod(N, f0, fp);
    SP = SPsetup(smod, [-1 2 -1]/sqrt(6), r_win{winix}(N-2));

    amp = 1;
    
    y = smod.s(amp, v);
    z = SP.F*y;
    
    % Prepare peak estimation array
    % We will use the fact that we know the true location.
    % Even the smallest peak will be identified. We do not try to
    % simulate distinguishability over noise.
    vest = ones(1,length(c_ip))*nan;

    for j=1:length(c_ip)
        [vests, ~, ~] = spectrum_peak_est(z,smod.vp,mph,c_ip(j));
        if ~isempty(vests)
            % Get closest peak
            [~,minix] = min(abs(vests - v));
            vest(1,j) = vests(minix);
        end
    end

    b(row,:) = mean(vest - v, 1, 'omitnan');

    fprintf(f, "%d & %s & %f & \n\t %f & %f & %f \\\\\n",...
        N, s_winname{winix}, v, b(row,:));
    
    fprintf('%d/%d in %f\n', row, R, toc(timer));
end

fclose(f);
