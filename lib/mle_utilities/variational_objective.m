function [f, gradf, hessf, t] = variational_objective(X, T, grad)
% Number of nodes in system.
nnodes = size(X, 2);

% Calculate the empirical mean parameters.
mu = mean_parameters(X);

% Dichotomous variable for partition function.
x = dec2bin(0:2^nnodes - 1) - '0'; % x(x == 0) = -1;

% Set up a symbolic matrix of real-valued random variables that correspond
% to the edge weights in the network. Set the lower triangular portion and
% the main diagonal to zero as we will not be needing these components.
omega = sym('x', nnodes);
omega = triu(sym(omega, 'real'));
omega(1:size(omega, 1) + 1:size(omega, 1)^2) = 0;

% keyboard
% This is also for structural knowledge.
% omega = omega .* (T ~= 0);


% We are now going to build up a symbolic function that represents that the
% objective function that we will try to optimize in order to recover the
% edge weight parameters.
f = 0;
partition = 0;

% The first component of the optimization function is just the product of
% the edge weights with their corresponding mean-value parameters. We loop
% through each vertex in the graph and calculate this symbolically.
for i = 1:nnodes
    for j = i:nnodes
        if i ~= j
            f = f + omega(i, j) * mu(i, j);
        end
    end
end

% The second component of objective function is the logarithm of the
% partition function. This means that we need to construct a symbolic
% representation of the partition function. We loop first through every
% configuration of the random vector for the Ising model.
for i = 1:size(x, 1)
    exponent = 0;
    % And then for each vertex, we calculate the pairwise product with
    % another vertex and accumulate.
    for j = 1:nnodes
        for k = 1:nnodes
            exponent = exponent + omega(j, k) * x(i, j) * x(i, k);
        end
    end

    partition = partition + exp(exponent);
end

% MATLAB solves minimization functions, so we convert our maximization
% function into a minimization by taking the negative.
f = -(f - log(partition));

% We construct a vector of the edge weight coefficients. Notice that we are
% only taking the upper triangular portion of the matrix (essentially).
t = omega(:);
t = t(t ~= 0);

% Calculate the Jacobian and the Hessian matrix of the symbolic
% representation of the optimization function.
gradf = jacobian(f, t).';

if grad
    hessf = [];
else
    hessf = jacobian(gradf, t);
end
