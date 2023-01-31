function dbplot(y,vp)
N = length(y);
plot(vp*(0:N-1)/N - vp/2, fftshift(20 .* log10(abs(y))), 'Color', [0 0 1]);
xlim([-vp/2 vp/2]);
xlabel('Velocity (m/s)')
ylim([min([ylim -20]) max([ylim 80])]);
ylabel('Amplitude (dB)')
grid on

%ax2 = axes('Position',gca().Position,'XAxisLocation','top','YAxisLocation','right','Color','none','YTick',[],'TickLength',[.00014 .0014]);
%xlabel(ax2,'Frequency (Hz)');

end
