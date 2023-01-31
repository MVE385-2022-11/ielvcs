% Frequency responses for MTI filters
clf

h_0 = [1];
h_1 = [-1 1]/sqrt(2);
h_2 = [-1 2 -1]/sqrt(6);
h_2_tilde = normalize([-1 16 -30 16 -1], 'norm', 2);  % Experiment
h_4 = normalize([1 -4 6 -4 1], 'norm', 2);  % Experiment

fp = 1;
fs = linspace(0,fp,1000);
clf
hold on
grid on

viridis_map = [
    253, 231, 37;
    94, 201, 98;
    33, 145, 140;
    59, 82, 139;
    68, 1, 84;
    ]/255;

cellfun(@(h,col) plot(fs, 20*log10(abs(freqz(h,1,fs,fp))), 'Color', col),...
    {h_0, h_1, h_2, h_4, h_2_tilde},...
    mat2cell(viridis_map,ones(1,size(viridis_map,1))).'...
    );

legend('h_0','h_1','h_2','h_4',[char([double('h') 771]) '_2'], 'Location', 'South');
ylabel('Amplitude (dB)');
xlabel('Doppler frequency (units of $$f_p$$)');
ylim([-50 10]);
