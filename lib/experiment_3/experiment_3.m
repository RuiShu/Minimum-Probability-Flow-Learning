function experiment_3(nsamples, node_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factorized graph estimation by MPF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings for MSE computation
ntrials = 50;

%% Settings for data generation
grid_shape = [1 node_size];      
file = ['experiment_3/experiment_3_n' int2str(nsamples) '_nodes' ...
        int2str(node_size) '.mat'];
fprintf('To be saved in: %s\n', file);

%% Generate weights
adj = ones(node_size) - eye(node_size);
J = ising_edge_generation(adj);

% burnin = 100*size(adj, 1);
% independent_steps = 10*size(adj, 1);
burnin = 1000;
independent_steps = 300;

%% Initialize data-recording variables
%% recorder is (ntrials, 6) matrix:
%% recorder(:,1) = nsamples
%% recorder(:,2) = node_size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf sq error
%% recorder(:,6) = mpf time
recorder = zeros(ntrials, 6);

fprintf('Data info: Nodes: %4d. Samples: %8d\n', ...
        size(adj, 1), nsamples);

for trial = 1:ntrials % can use parfor loop for parallelization
    Xall_bin = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps);

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
    recorder(trial, 2) = node_size;
    recorder(trial, 3) = mle_sq_err;
    recorder(trial, 4) = mle_t;
    recorder(trial, 5) = mpf_sq_err;
    recorder(trial, 6) = mpf_t;

    %% Display the values
    fprintf(['n: %5d, nodes: %3d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
             '%10.4f, mpf time: %10.4f\n'], ...
            nsamples, node_size, mle_sq_err, mle_t, mpf_sq_err, mpf_t);
end

%% Results
results = mean(recorder, 1);

fprintf('***Final Analysis***\n');
fprintf(['n: %5d, nodes: %3d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
         '%10.4f, mpf time: %10.4f\n'], ...
        results(1), results(2), results(3), results(4), results(5), ...
        results(6));

%% Save results
save(file, 'recorder', 'results');
