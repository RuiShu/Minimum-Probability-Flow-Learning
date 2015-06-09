function experiment_1(nsamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factorized graph estimation by MPF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings for MSE computation
ntrials = 30;

%% Settings for data generation
grid_shape = [4 4];                     % Do not change
file = ['experiment_1/experiment_1_n' int2str(nsamples) '.mat'];
fprintf('To be saved in: %s\n', file);

%% Generate weights
adj = adj_mat_generation(grid_shape); % adjacency matrix
J = ising_edge_generation(adj);

% burnin = 100*size(adj, 1);
% independent_steps = 10*size(adj, 1);
burnin = 1000;
independent_steps = 300;

%% Initialize data-recording variables
%% recorder is (ntrials, 5) matrix:
%% recorder(:,1) = nsamples
%% recorder(:,2) = mle sq error
%% recorder(:,3) = mle time
%% recorder(:,4) = mpf sq error
%% recorder(:,5) = mpf time
recorder = zeros(ntrials, 5);

fprintf('Data info: Nodes: %4d. Samples: %8d\n', ...
        size(adj, 1), nsamples);

for trial = 1:ntrials % can use parfor loop for parallelization
    Xall_bin = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps); % for MPF

    fprintf('Running analysis\r');
    %% Maximum Likelihood Estimation
    mle_t = tic();
    J_estimate = mle_estimate(adj, Xall_bin);
    mle_t = toc(mle_t);
    mle_sq_err = norm(J(:) - J_estimate(:))^2;

    %% Minimum Probability Flow
    mpf_t = tic();
    J_estimate = mpf_estimate(adj, Xall_bin);
    mpf_t = toc(mpf_t);
    mpf_sq_err = norm(J(:) - J_estimate(:))^2;

    %% Record the values
    recorder(trial, 1) = nsamples;
    recorder(trial, 2) = mle_sq_err;
    recorder(trial, 3) = mle_t;
    recorder(trial, 4) = mpf_sq_err;
    recorder(trial, 5) = mpf_t;

    %% Display the values
    fprintf(['n: %5d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
             '%10.4f, mpf time: %10.4f\n'], ...
            nsamples, mle_sq_err, mle_t, mpf_sq_err, mpf_t);
end

%% Results
results = mean(recorder, 1);

fprintf('***Final Analysis***\n');
fprintf(['n: %5d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
         '%10.4f, mpf time: %10.4f\n'], ...
        results(1), results(2), results(3), results(4), results(5));

%% Save results
save(file, 'recorder', 'results');
