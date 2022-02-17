function hazard = harzard_cle(x, yT, sys, dt, c)
%--------------------------------------------
% h(x_t|y_T) = a_j(x_t)*pcle(y_T | X_t = x_t + nu_j)/pcle(y_T | X_t = x_t).

% where pcle(y_T | X_t = x_t) = Normalpdf(y_T, mu, Sigma).

% with mu = Proj*(x_t + nu*a(x_t)*dt), 
% where Proj is the projection mapping, Y = Proj(X, Y)
% i.e. The expected V(T) for a tau-leaping process
% with dt = tau = T- t, h(s) = a(x_t) for t<s<T).


% with Sigma = Proj*nu*H(x_t)*nu'*Proj'*dt
% where H(x_t) = diag(a(x_t)).
% ---------------------------------------------

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n2 = length(yT);
Proj = [zeros(n2, n-n2), eye(n2)]; 

lambda = feval(sys, 'prop', x, c);
hazard = zeros(m, 1);

mu = Proj*(x + nu*lambda*dt);
H = diag(lambda);
Sigma = Proj*nu*H*nu'*Proj'*dt + 5*eye(2);
%Sigma = Proj*nu*H*nu'*Proj'*dt;

pcle_xt = mvnpdf(yT, mu, Sigma);
%mvnpdf: multi-variate normal pdf

for j=1:m
    x_next = x + nu(:,j);
    lambda_next = feval(sys, 'prop', x_next, c);
    mu = Proj*(x_next + nu*lambda_next*dt);
    H = diag(lambda_next);
    Sigma = Proj*nu*H*nu'*Proj'*dt;
    pcle_x_j = mvnpdf(yT, mu, Sigma);
    hazard(j) = lambda(j)*pcle_x_j/pcle_xt;
end
end