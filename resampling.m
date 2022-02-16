function [V_new, w_new, k_new] = resampling(V, w, k)
    wsum = sum(w);
    w = w/wsum;
    [n_unobs, Ns] = size(V);
    [m, Ns] = size(k);
    V_new = zeros(n_unobs, Ns);
    w_new = zeros(1, Ns);
    k_new = zeros(m, Ns);
    
    offs = offsprings(w);
    %sort_offs = sort(offs, 'descend');
    i=1; ind=1;
    for i=1:Ns
        w_new(i)=1;
        for l=1:offs(i)
            V_new(:,ind) = V(:,i);
            k_new(:,ind) = k(:,i);
            ind = ind+1;
        end
    end
end