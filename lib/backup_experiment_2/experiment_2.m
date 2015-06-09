function experiment_2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factorized graph estimation by MPF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings for MSE computation
ntrials = 10;

%% Settings for data generation
grid_shape = [15 15];
n = [1, 3, 5, 15];

%% Check that grid_shape is square
if grid_shape(1) ~= grid_shape(2)
    fprintf('Not a square\n');
    return
end

%% Check that all values in n are factors of grid_shape dim
if any(bsxfun(@mod, grid_shape(1), n) > 0)
    fprintf('Not a value n_size partitioning\n');
    return
end

nsamples = 1000;
adj = adj_mat_generation(grid_shape); % adjacency matrix

% burnin = 100*size(adj, 1);
% independent_steps = 10*size(adj, 1);
burnin = 1000;
independent_steps = 300;

%% Initialize data-recording variables
sq_err_MPF = zeros(ntrials, length(n));
time_MPF = zeros(ntrials, length(n));
J = ising_edge_generation(adj);

for trial = 1:ntrials % can use parfor loop for parallelization
    Xall_bin = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps); % for MPF

    fprintf('Data generated. Nodes: %4d. Samples: %8d:\n', ...
            size(Xall_bin, 1), size(Xall_bin, 2));

    %% Minimum Probability Flow
    for n_size = n
        node_set = get_node_set(grid_shape, n_size, adj);
        J_estimate = zeros(size(adj));

        t = tic();
        for i = 1:size(node_set, 1)
            inner_nodes = node_set{i, 1};
            nodes = node_set{i, 2};

            J_partition = mpf_multiedge_estimate(adj, Xall_bin, ...
                                                 inner_nodes, nodes);

            %% Print partial error
            est = J_partition;
            true = J(nodes, nodes);
            fprintf(': %.4f :\r', norm(true(:) - est(:))^2);

            %% Store in J_estimate
            J_estimate(nodes, nodes) = J_partition;

        end

        t = toc(t);
        sq_err = norm(J(:) - J_estimate(:))^2;

        %% Record the values
        time_MPF(trial, find(n == n_size)) = t;
        sq_err_MPF(trial, find(n == n_size)) = sq_err;

        %% Display the values
        fprintf('Square error: %2.4f, Time: %2.4f\n', sq_err, t);
    end
end

%% Results
fprintf('MPF Average run time:      ');
disp(mean(time_MPF))

fprintf('MPF MSE:                   ');
disp(mean(sq_err_MPF));
