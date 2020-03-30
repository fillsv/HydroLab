function ind = findnear(data, value)
    ind = find(min(abs(data-value))==abs(data-value));
    ind = ind(1);
    