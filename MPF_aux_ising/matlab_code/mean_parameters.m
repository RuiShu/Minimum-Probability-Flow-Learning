function mu = mean_parameters(X)
nnodes = size(X, 2);
mu = zeros(nnodes);
for i = 1:nnodes
    for j = i+1:nnodes
        mu(i,j) = mean(X(:, i) .* X(:, j));
    end 
end
