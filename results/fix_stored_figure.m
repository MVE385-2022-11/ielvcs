if exist('sim_grid', 'var') ~= 1
    error('No data loaded');
end
fprintf('Current row: %d\n', row);

f = gcf;
t = f.Children(1);

width = 5.93;%3.3;     % Width in inches
height = 3;%3.3;    % Height in inches
f.PaperSize = [width height];
f.PaperPosition = [0 0 width height];

%%
t.Children(6).Title.String = sprintf("SP: b = %.3f ± %.3f (%.2f)",...
    b(row,1), b_err(row,1), b_sucfrac(row,1));
t.Children(5).Title.String = sprintf("MF: b = %.3f ± %.3f (%.2f)",...
    b(row,2), b_err(row,2), b_sucfrac(row,2));
t.Children(4).Title.String = sprintf("MF (os): b = %.3f ± %.3f (%.2f)",...
    b(row,3), b_err(row,3), b_sucfrac(row,3));
t.Children(2).Title.String = sprintf("LASSO: \x03bb = %g, b = %.3f ± %.3f (%.2f)",...
    lambdas(1,winix), b(row,4), b_err(row,4), b_sucfrac(row,4));
t.Children(1).Title.String = sprintf("LASSO (os): \x03bb = %g, b = %.3f ± %.3f (%.2f)",...
    lambdas(2,winix), b(row,5), b_err(row,5), b_sucfrac(row,5));
for i=1:6
    t.Children(i).Title.FontSize = 6.5;
    t.Children(i).XLim = [-25 115];
end

%% one_target
winix = sim_grid(row, 2);
snr = r_SNR(sim_grid(row, 3));
v = sim_grid(row, 4);

set(f, 'CurrentAxes', t.Children(3))
cla
axis ij
xlim([0 1]);
set(gca,'Visible','off');
text(0,0.5,sprintf("K = %d\nN = %d\nN' = [%s]\nv = %.3f m/s\nSNR = %.1f dB\nwindow = %s",...
    K, N, sprintf("%d, %d", Nprime1, Nprime2), v, snr, s_winname{winix}), 'FontSize', 9);

%% separation
gr = num2cell(sim_grid(row,:));
[N, winix, snrix, v1, dv] = deal(gr{:});
snr = r_SNR(snrix);
v2 = v1 + dv;

t.Children(6).Title.String = sprintf("SP: %% = %.3f ± %.3f",...
    sucfrac(row,1), sucfrac_err(row,1));
t.Children(5).Title.String = sprintf("MF: %% = %.3f ± %.3f",...
    sucfrac(row,2), sucfrac_err(row,2));
t.Children(4).Title.String = sprintf("MF (os): %% = %.3f ± %.3f",...
    sucfrac(row,3), sucfrac_err(row,3));
t.Children(2).Title.String = sprintf("LASSO: \x03bb = %g, %% = %.3f ± %.3f",...
    lambdas(1,winix), sucfrac(row,4), sucfrac_err(row,4));
t.Children(1).Title.String = sprintf("LASSO (os): \x03bb = %g, %% = %.3f ± %.3f",...
    lambdas(2,winix), sucfrac(row,5), sucfrac_err(row,5));
for i=1:6
    t.Children(i).Title.FontSize = 6.5;
    t.Children(i).XLim = [-5 20];
    t.Children(i).YLim = [-10 60];
end

set(f, 'CurrentAxes', t.Children(3))
cla
axis ij
xlim([0 1]);
ylim([0 1]);
set(gca,'Visible','off');
text(0,0.5,sprintf("K = %d\nN = %d\nv1 = %.3f\nv2 = %.3f\ndv = %.3f\nSNR = %.3f\nwindow = %s",...
            K, N, v1, v2, dv, snr, s_winname{winix}), 'FontSize', 9);
