function mu = mean_parameters(X)

disp(X)

mu = zeros(size(X, 2));

for i = 1:size(X, 2)
    for j = (i+1):size(X, 2)
        mu(i, j) = mean(X(:, i).*X(:, j));
        mu(j, i) = mu(i, j);
    end
end

disp(mu)

