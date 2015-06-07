function experiment_one()

%% Settings for MSE computation
ntrials = 1;

%% Settings for data generation
grid_shape = [30 30];
n = 3:28;
nsamples = 10000;
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

    fprintf(['Finished generating data Xall. size(Xall) = [# nodes, ' ...
             '# samples]:\n']);
    disp(size(Xall_bin));
    
    %% Minimum Probability Flow
    for n_size = n
        J_estimate = zeros(n_size^2);
        [inner_nodes, nodes] = get_all_nodes(grid_shape, n_size, adj);

        t = tic();
        [rui_inner_idx, J_estimate] = mpf_centered_square_estimate(adj, Xall_bin, ...
                                                  inner_nodes, nodes);
        t = toc(t);

        est = J_estimate;
        true = J(nodes, nodes);

        true = true(rui_inner_idx, rui_inner_idx);
        est = est(rui_inner_idx, rui_inner_idx);

        disp(norm(true - est))

        fprintf('Done\n');
%         keyboard;

%         time_MPF(trial) = t;
%         fprintf('Square error: ');
%         sq_err_MPF(trial) = sum((J(:) - J_estimate(:)).^2);
%         sq_err_Null(trial) = sum(J(:).^2);
%         disp(sq_err_MPF(trial));
    end
end

%% Results
fprintf('MPF Average run time:      ');
disp(mean(time_MPF))

fprintf('MPF MSE:                   ');
disp(mean(sq_err_MPF));
