function main()

%% Settings for MSE computation
ntrials = 1;
sq_err_MPF = zeros(ntrials, 1);
sq_err_Null = zeros(ntrials, 1);
time_MPF = zeros(ntrials, 1);

%% Settings for data generation
grid_shape = [2, 3];
nsamples = 10000;
adj = adj_mat_generation(grid_shape); % adjacency matrix

J = ising_edge_generation(adj);

for trial = 1:ntrials % can use parfor loop for parallelization
    Xall = ising_data_generation(J, grid_shape(1)*grid_shape(2), ...
                                 nsamples);
    fprintf(['Finished generating data Xall. size(Xall) = [# nodes, ' ...
             '# samples]:\n']);
    disp(size(Xall));

    %% Minimum Probability Flow
    fprintf('MPF Estimation:\n');
    J_estimate = zeros(size(adj));

    t = tic();
    for j = 1:size(J, 1)
        for k = (j+1):size(J, 2)
            J_estimate(j, k) = ising_edge_estimation(Xall, adj, j, k);
            J_estimate(k, j) = J_estimate(j, k);
        end
        fprintf('%.3f complete. Trial %3d out of %3d\n', ...
                j/size(J, 1), trial, ntrials);
    end
    t = toc(t);
    time_MPF(trial) = t;
    sq_err_MPF(trial) = sum((J(:) - J_estimate(:)).^2);
    sq_err_Null(trial) = sum(J(:).^2);

    fprintf('Square error: ');
    disp(sum((J(:) - J_estimate(:)).^2));

    %% Maximum Likelihood Estimation
    fprintf('MLE Estimation:\n');
    [theta_est, e] = variational_estimate(J, Xall', ...
                                          1, 1:size(J, 1));
    fprintf('Square error: ');
    disp(e);
end

%% Results
fprintf('Average run time:    ');
disp(mean(time_MPF))

fprintf('MPF MSE:             ');
disp(mean(sq_err_MPF));

fprintf('Null hypothesis MSE: ');
disp(mean(sq_err_Null));
