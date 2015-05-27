function p = unnormalized_probability(x, T)
d = size(T, 1);
p = 0;
for i = 1:d
    for j = i+1:d
        p = p + T(i, j) * x(i) * x(j);
    end
end
p = exp(p);