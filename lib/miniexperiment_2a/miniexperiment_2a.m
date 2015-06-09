function miniexperiment_2a()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double patch experiment for MPF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings for data generation
grid_shape = [10 10];
n = 3:28;
nsamples = 10000;
adj = adj_mat_generation(grid_shape); % adjacency matrix

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
    %% Patch analysis for 10 by 10
    %% Patch one:
    fprintf('Patch One:\n');
    inner_nodes = [23 24 33 34]';
    nodes = [23 24 33 34 13 14 22 25 32 35 43 44]';
    [rui_inner_idx, J_estimate] = mpf_multiedge_estimate(adj, Xall_bin, ...
                                                      inner_nodes, nodes);
    est = J_estimate;
    true = J(nodes, nodes);
    disp(norm(true - est))

    %% Patch two:
    fprintf('Patch Two:\n');
    inner_nodes = [68 69 78 79]';
    nodes = [68 69 78 79 58 59 67 70 77 80 88 89]';
    [rui_inner_idx, J_estimate] = mpf_multiedge_estimate(adj, Xall_bin, ...
                                                      inner_nodes, nodes);
    est = J_estimate;
    true = J(nodes, nodes);
    disp(norm(true - est))

    %% Patch Three:
    fprintf('Patch Three:\n');
    inner_nodes = [23 24 33 34 68 69 78 79]';
    nodes = [23 24 33 34 13 14 22 25 32 35 43 44 68 69 78 79 58 59 67 70 77 80 88 89]';
    [rui_inner_idx, J_estimate] = mpf_multiedge_estimate(adj, Xall_bin, ...
                                                      inner_nodes, nodes);
    est = J_estimate;
    true = J(nodes, nodes);
    disp(norm(true - est))

    fprintf('Done\n');

end
