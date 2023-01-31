function [beta, n] = LASSO_ISTA(SPy, beta0, lambda, G, step, tolerance)

% objective function
obj = @(beta) norm(G*beta - SPy)^2 + lambda*norm(beta,1);

% soft thresholding function
eta = @(z, gamma) max(0, abs(z) - gamma).*sign(z);

t = step;           % stepsize in (0, 1/(2*g2_norm));
tol = tolerance;    %7e2 * g2_norm;   

GtG = G'*G;
GtS = G'*SPy;

n = 1;
beta_last = beta0;  % Init beta
beta = eta(beta_last - t*(GtG*beta_last - GtS), t*lambda);
while abs(obj(beta)-obj(beta_last)) > tol
    beta_last = beta;
    beta = eta(beta_last - t*(GtG*beta_last - GtS), t*lambda);
    n = n+1;
    if mod(n,1000) == 0
        fprintf(' ...%d @ delta %f\n', n, abs(obj(beta)-obj(beta_last)))
    end
end