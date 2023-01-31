N = 1024+2;
f0 = 3e9; 
fp = 1e4;


smod = sigmod(N, f0, fp);
win{1} = @hamming;
win{2} = @hann;
win{3} = @welch;
win_names = ["Hamming" "Hann" "Welch"];


%% Show equalizer frequency response
clf
hold on
xlabel('Velocity (m/s)');
ylabel('Amplitude (dB)');
yline(0, 'k:');

for j = 1:length(win)
    SP = SPsetup(smod,[-1 2 -1]/sqrt(6), win{j}(N-2));
    e = vecnorm(SP.C,2,2);

    plot(smod.vp*(0:SP.M-1)/SP.M - smod.vp/2, fftshift(20*log10(abs(1./e))));
end
legend(["" win_names]);
xlim([-30 30]);
ylim([0 80]);


