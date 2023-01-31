function [beta, n] = LASSO_FISTA(SPy, beta0, lambda, G, tol)
% Beck, Teboulle : A Fast Iterative Shrinkage-Thresholding Algorithm for
% Linear Inverse Problems

%scaling?
%scale = norm(SPy,2)^2;
%lambda = lambda.*norm(SPy,1);
%tol = tol.*norm(SPy,2)^2;

% objective function
obj = @(beta) norm(G*beta - SPy)^2 + lambda*norm(beta,1);

% soft thresholding function
eta = @(z, gamma) max(0, abs(z) - gamma).*sign(z);

% gradient of smooth term
GtG = G'*G;
GtS = G'*SPy;
grad = @(beta) 2*(GtG*beta - GtS); %scale?
L = 2*max(eig(GtG)); %scale?
Linv = 1/L;

n = 1;
beta_y = beta0;  % Init beta
beta_x = beta0;
beta_x_last = beta0;
t = 1;

while n < 10 || (abs(obj(beta_x)-obj(beta_x_last)) > tol)
    beta_x_last = beta_x;
    beta_x = eta(beta_y - Linv*(GtG*beta_y - GtS), Linv*lambda);

    t_last = t;
    t = (1 + sqrt(1 + 4*t_last^2))/2;

    beta_y = beta_x + (t_last - 1)/t * (beta_x - beta_x_last);

    n = n+1;
    if mod(n,1000) == 0
        fprintf(' ...%d delta %.1f\n', n, abs(obj(beta_x)-obj(beta_x_last)))
    end
end

beta = beta_x;

end