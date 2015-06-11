function node_set = get_chimera_node_set(n_size, adj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keywords:                                                                                
% n_size     -- A numeric describing number of inner nodes (batches of chimeras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:(size(adj, 1)/n_size)
    [inner_nodes, nodes] = get_nodes_by_index((i-1)*n_size + 1, n_size, adj);
    node_set{i, 1} = inner_nodes;
    node_set{i, 2} = nodes;
end
