function miniexperiment_2b(n_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factorized graph estimation by MPF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings for MSE computation
ntrials = 10;

%% Settings for data generation
nsamples = 1000;
grid_shape = [30 30];
save_file = ['miniexperiment_2b/miniexperiment_2b_factorization' int2str(n_size) '.mat'];
fprintf('To be saved in: %s\n', save_file);

%% Generate adj
adj = adj_mat_generation(grid_shape); % adjacency matrix

%% Initialize data-recording variables
%% recorder is (ntrials, 6) matrix:
%% recorder(:,1) = sample size
%% recorder(:,2) = factorization size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf sq error
%% recorder(:,6) = mpf time
recorder = zeros(ntrials, 6);


dataset = load('data/dataset1.mat');
Xall_bin = datasample(dataset.Xall, nsamples, 2, 'Replace', false);

fprintf('Data info: Nodes: %4d. Samples: %8d\n', ...
        size(Xall_bin, 1), size(Xall_bin, 2));

fprintf('Generating node_set\r');
node_set = get_node_set(grid_shape, n_size, adj);

for trial = 1:ntrials % can use parfor loop for parallelization
    %% Get dataset
    load_file = ['data/dataset' int2str(trial) '.mat'];
    fprintf('Loading dataset, %s\r', load_file);

    dataset = load(load_file);
    Xall_bin = datasample(dataset.Xall, nsamples, 2, 'Replace', false);
    J = dataset.J;

    %% Maximum Likelihood Estimation
    if n_size > 2
        mle_sq_err = nan;
        mle_t = nan;
    else 
        J_estimate = zeros(size(adj));

        mle_t = tic();
        for i = 1:size(node_set, 1)
            inner_nodes = node_set{i, 1};
            nodes = node_set{i, 2};

            J_partition = mle_multiedge_estimate(adj, Xall_bin, ...
                                                 inner_nodes, nodes);

            %% Print partial error
            est = J_partition;
            true = J(nodes, nodes);
            fprintf(': %.4f : mle, %3d\r', norm(true(:) - est(:))^2, i);

            %% Store in J_estimate
            J_estimate(nodes, nodes) = J_partition;

        end
        mle_t = toc(mle_t);
        mle_sq_err = norm(J(:) - J_estimate(:))^2;
    end

    %% Minimum Probability Flow
    J_estimate = zeros(size(adj));

    mpf_t = tic();
    for i = 1:size(node_set, 1)
        inner_nodes = node_set{i, 1};
        nodes = node_set{i, 2};

        J_partition = mpf_multiedge_estimate(adj, Xall_bin, ...
                                             inner_nodes, nodes);

        %% Print partial error
        est = J_partition;
        true = J(nodes, nodes);
        fprintf(': %.4f : mpf, %3d\r', norm(true(:) - est(:))^2, i);

        %% Store in J_estimate
        J_estimate(nodes, nodes) = J_partition;

    end
    mpf_t = toc(mpf_t);
    mpf_sq_err = norm(J(:) - J_estimate(:))^2;

    %% Record the values
    recorder(trial, 1) = nsamples;
    recorder(trial, 2) = n_size;
    recorder(trial, 3) = mle_sq_err;
    recorder(trial, 4) = mle_t;
    recorder(trial, 5) = mpf_sq_err;
    recorder(trial, 6) = mpf_t;

    %% Display the values
    fprintf(['n: %5d, n_size: %3d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
             '%10.4f, mpf time: %10.4f\n'], ...
            size(Xall_bin, 2), n_size, mle_sq_err, mle_t, mpf_sq_err, mpf_t);

end

%% Results
results = mean(recorder, 1);

fprintf('***Final Analysis***\n');
fprintf(['n: %5d, n_size: %3d, mle sq_err: %10.4f, mle time: %10.4f, mpf sq_err: ' ...
         '%10.4f, mpf time: %10.4f\n'], ...
        results(1), results(2), results(3), results(4), results(5), ...
        results(6));

%% Save results
save(save_file, 'recorder', 'results');
