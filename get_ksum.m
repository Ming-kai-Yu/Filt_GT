function ksum = get_ksum(w)

w = w/sum(w);
wsorted = sort(w, 'descend');
ksum = cumsum(wsorted);
end