function freq = get_hist(data, weight, bin)
% Get histcount of weighted data

num_bin = length(bin);
freq = zeros(1, num_bin);
weight = weight/sum(weight);
for i=1:num_bin
    ind = (data == bin(i));
    freq(i) = sum(weight(ind));
end

freq = freq/sum(freq);
end