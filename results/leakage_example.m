% Figure that illustrates leakage in DFT

Fs = 100; % sampling frequency
f0 = 5; % signal frequency
T = .75; % time duration 
t = 0:(1/Fs):T;  %# the time samples
t_fine = 0:(1/(100*Fs)):T;
N = length(t);
f = (-Fs/2):(Fs/N):(Fs/2);  % the frequency bin axis
f = f(1:end-1);

cos_t = cos(2*pi*f0*t);  % generate the sampled signal
cos_t_fine = cos(2*pi*f0*t_fine); 

%%
clf
stem(t, cos_t, 'xb');
hold on
plot(t_fine, cos_t_fine, 'Color', [0 0 1 0.15]);
xlim([0, 1]);

%%
clf
y_1 = fftshift(abs(fft(cos_t)));
stem(f, y_1, 'xb-', 'MarkerSize', 8);
xlim([-3 3]*f0);
ylim([0 60]);

%%
clf
y_1 = fftshift(abs(fft(hamming(length(cos_t)).' .* cos_t)));
stem(f, y_1, 'xb-', 'MarkerSize', 8);
xlim([-3 3]*f0);
ylim([0 60]);


