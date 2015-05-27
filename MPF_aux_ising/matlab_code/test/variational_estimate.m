function [theta_est, e] = variational_estimate(T, nsamples, niters, ...
                                               grad, indices)
% Load in the mean value parameters, which were either analytical or
% empirical. We will want to compare how using the analytical variety
% compares to the empirical results.
X = sample_ising(T, 10000, 900, 90)';

fprintf('finished generating samples\n');

T = T(indices, indices); X = X(:, indices);

% Reformat the edge weight matrix.
theta = cell(length(2:size(T, 1)), 1);
for i = 2:size(T, 1)
    theta{i - 1} = T(1:i-1, i);
end

% This yields a theta vector that can be directly compared to the vector
% that we generate through variational optimization.
theta = vertcat(theta{:});



% Create the variational objective function.
[f, gradf, hessf, t] = variational_objective(X, T, grad);

% Create a symbolic function that can be optimized with respect to the
% vector of edge weights.

if grad
    fh = matlabFunction(gradf, 'vars', {t});
else
    fh = matlabFunction(f, gradf, hessf, 'vars', {t});
end


% Solve the variational optimization problem that recovers the original
% edge weight parameters.

if grad
    % Initialize the starting point for gradient descent.
    theta_est = zeros(size(t));
    theta_temp = inf(size(t));
    
    % Set the learning rate and the extent to which the learning rate
    % decreases after each iteration.
    lr = 1;alpha = 0.99;
    iter = 0;
    while norm(theta_est - theta_temp) > 1e-16
        if mod(iter, 500) == 0
            fprintf('iteration: %d\n', iter);
        end
        iter = iter + 1;
        theta_temp = theta_est;        
        deriv = fh(theta_est);        
        theta_est = theta_est - lr * deriv;
        lr = lr * alpha;
    end
else
    % Set options for the variational optimization problem. The most
    % important properties are the maximum number of allowed function
    % evaluations and the maximum number of allowed iterations of the
    % algorithm.
    options = optimoptions('fminunc', 'GradObj', 'on', 'Hessian', 'on', ...
        'Display', 'off', 'TolFun', 1e-16, 'MaxIter', 2000, ...
        'MaxFunEvals', 2000);
    theta_est = fminunc(fh, zeros(nchoosek(size(T, 1), 2), 1), options);
end

% Construct a trimmed-down version of theta.
theta_trim = theta(theta ~= 0);
theta_est_trim = theta_est(theta ~= 0);

% fprintf('True edge weight parameter:\n');
% disp(theta_trim);
% 
% fprintf('Estimated edge weight parameter:\n');
% disp(theta_est_trim);

% Calculate the error.
e = (theta_trim - theta_est_trim);e = e'*e;

disp([theta_trim, theta_est_trim]);
disp(e);

fprintf('Done');