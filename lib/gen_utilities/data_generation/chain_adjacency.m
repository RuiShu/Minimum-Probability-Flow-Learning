function A = chain_adjacency(n_nodes)

A = zeros(n_nodes);A(2:n_nodes+1:end) = 1;A = A + A';
