function stability = dbplot_ALL(Y,vp)
Y = fftshift(20 .* log10(abs(Y)), 1);
N = size(Y,1);
K = size(Y,2);
vs = (vp*(0:N-1)/N - vp/2).';

inc = any(Y > -100, 2);
stability = mean(Y > -100, 2);

hold on
%stem(vs, Y, 'Color', [0 0 0 0.15]);
%plot(vs, mean(Y,2), 'k');
for k=1:K
    scatter(vs, Y(:,k), 3.5, [stability - (stability == 1) zeros(N,1) stability], 'AlphaData', 0.5);
end

% Stability text
%unst = stability < 1 & inc;
%text(vs(unst), max(Y(unst,:),[],2) + 10, sprintfc('%.2f',stability(unst)), 'Rotation', 90);

xlim([-vp/2 vp/2]);
xlabel('Velocity (m/s)')
ylim([min([ylim -20]) max([ylim 80])]);
ylabel('Amplitude (dB)')
grid on

end
