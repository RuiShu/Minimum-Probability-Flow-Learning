function generate_chimera_data()

%% Grid info for ising data
n_rows = 5;
n_columns = 5;
n_nodes = 4;
nsamples = 10000;

adj = chimera_adjacency(n_rows, n_columns, n_nodes);

%% Generate chimera data
for i = 1:100
    file = ['data/chimera_dataset' int2str(i) '.mat'];
    fprintf('%s\n', file);

    J = ising_edge_generation(adj);

    burnin = 1000;
    independent_steps = 300;

    Xall = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps);

    %% Save results
    save(file, 'J', 'Xall');

end
