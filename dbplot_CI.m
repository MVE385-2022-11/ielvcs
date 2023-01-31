function dbplot_CI(Y,vp)
Y = 20 .* log10(abs(Y));
N = size(Y,1);

y_1 = Y(:,1);
Q = quantile(Y.', [0.025 0.5 0.975]).';

y_left = Q(:,1);
y_median = Q(:,2);
y_right = Q(:,3);

vs = (vp*(0:N-1)/N - vp/2).';

hold on
fill([vs.' fliplr(vs.')], [fftshift(y_left).' fliplr(fftshift(y_right).')],...
    'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
plot(vs, fftshift(y_1), 'Color', [0 0 1 0.15]);
plot(vs, fftshift(y_median), 'k');

xlim([-vp/2 vp/2]);
xlabel('Velocity (m/s)')
ylim([min([ylim -20]) max([ylim 80])]);
ylabel('Amplitude (dB)')
grid on

end
