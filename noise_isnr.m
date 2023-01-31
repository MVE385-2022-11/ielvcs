function sigma=noise_isnr(snr, N, amp)
% NOISE_ISNR  Computes the correct noise variance to achieve a given
%             coherently integrated SNR in a single target setting.
%   sigma = NOISE_ISNR(snr, N, amp) computes the variance for
%     SNR, number of pulses, and target amplitude.

% SNR = 10*log10((s'*s)^2 / s'*s sigma^2) = 10*log10(s'*s / sigma^2)
%     = 10*log10(N amp^2/sigma^2) = 10*log10(N) + 20*log10(amp) - 20*log10(sigma)
% Hence
% sigma = 10^[ -SNR/20 + log10(amp) + log10(N)/2 ]
%       = amp*sqrt(N)*10^[-SNR/20]

sigma = amp * sqrt(N) * 10^(-snr/20);
end