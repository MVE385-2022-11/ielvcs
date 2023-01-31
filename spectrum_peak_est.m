function [vests, vlocs, locs] = spectrum_peak_est(z, vp, mph, k)
N = size(z,1);
vs = vp*(0:N-1)/N - vp/2;
k = (k-1)/2; % total pts to pts on each side

z = fftshift(20*log10(abs(z)));

[~,locs] = findpeaks(z, 'MinPeakHeight', mph);

% Peak locations by index
vlocs = vs(locs);

% Peak locations by quadratic fit using neighbours
vests = zeros(size(vlocs));
for i = 1:length(vlocs)
    % index cyclically if at edges
    xs = vs(1 + mod(locs(i) + (-k:k) - 1, N));
    ys = z(1 + mod(locs(i) + (-k:k) - 1, N));
    %xs = vs(locs(i) + (-k:k));
    %ys = z(locs(i) + (-k:k));
    ps = polyfit(xs,ys,2);
    vests(i) = -ps(2)/(2*ps(1));
end

end
