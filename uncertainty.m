%% Verification
% This script currently requires the results after running 
% the script 'run_four_species_full_lattice.m'.

%% L2 error

% x1_hat_dat(i,j) stores the expectation of X1 given X3 = i-1 and X4 = j-1
x1_hat_dat = zeros(base, base);
for i = 1:base
    for j = 1:base
        x1_hat_dat(i,j)=e_x1_given_x3x4(i-1, j-1, p_final, base);
    end
end

l2 = 0;
for i = 1:length(p_final)
    x = ind2state(i, base);
    l2 = l2 + (x(1) - x1_hat_dat(x(3)+1, x(4)+1))^2*p(i);
end




%% local functions

function x1_hat = e_x1_given_x3x4(x3, x4, p, base)

for i = 1:length(p)
    x = ind2state(i, base);
    p34 = 0;
    x1_hat = 0;
    if x(3) == x3 && x(4) == x4
        x1_hat = x1_hat + x(1)*p(i);
        p34 = p34 + p(i);
    end
    
    if p34 ~= 0
        x1_hat = x1_hat/p34;
    end
end

end %func
