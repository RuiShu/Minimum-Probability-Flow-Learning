function generate_ising_data()

%% Grid info for ising data
grid_shape = [30 30];
nsamples = 10000;

%% Generate 
adj = adj_mat_generation(grid_shape); % adjacency matrix

for i = 1:1000
    file = ['data/dataset' int2str(i) '.mat'];
    fprintf('%s\n', file);

    J = ising_edge_generation(adj);

    burnin = 1000;
    independent_steps = 300;

    Xall = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps);

    %% Save results
    save(file, 'J', 'Xall');

end
