L = 128;
M = 512;
xs = 0:L-1;
fs = (0:M-1) - M/2;

names = {'Hamming', 'Hann', 'Welch'};
windows = [hamming(L) hann(L) welch(L)];

t = tiledlayout(2, size(windows,2));
t.TileSpacing = 'compact';
t.Padding = 'compact';

for j = 1:size(windows,2)
    nexttile;
    plot(xs, windows(:,j), 'b');
    xlim([0 L-1]);
    ylim([0 1.1]);
    title(names{j});
    if j == 1
        ylabel('Time domain');
    end
    grid on;
end

for j = 1:size(windows,2)
    nexttile;
    fy = 20*log10(abs(fft(windows(:,j), M)));
    fy = fy - max(fy);
    stem(fs, fftshift(fy), 'Color', [1 0.5 0], 'Marker','none', 'BaseValue', -200);
    xlim([-100 100]); %xlim([-M/2 M/2]);
    ylim([-100 10]);
    if j == 1
        ylabel('Frequency domain');
    end
    grid on;
end
