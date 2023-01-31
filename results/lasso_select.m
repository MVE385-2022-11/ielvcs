clear

N = 130;
Nprime = 1024;
win = @hann;
amp = [1];
vs = [5];
sigma = noise_isnr(20,N,amp(1));

% General params
f0 = 3e9; 
fp = 1e4;
lambda_ref = 30; % only a reference

tol = 1e-4; % FISTA tolerance

smod = sigmod(N, f0, fp);
SP = SPsetup(smod, [-1 2 -1]/sqrt(6), win(N-2));

s = smod.s(amp, vs);
if sigma == 0
    n = zeros(smod.N,1);
else
    n = smod.randcn(sigma, smod.N, 1);
end
y = s + n;

[lambdas, paths] = lambda_paths(y, smod, SP, Nprime, [1 1000], tol);
pause(0.5);
lambda_color(vs, smod.vp/smod.N);

% Reference line
xline(lambda_ref,':',Parent=gca);
