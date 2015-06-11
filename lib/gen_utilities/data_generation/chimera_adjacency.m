function adj = chimera_adjacency(n_rows, n_columns, n_nodes)
I_rows = eye(n_rows);I_columns = eye(n_columns);I_nodes = eye(n_nodes);
O_nodes = ones(n_nodes);
L_rows = chain_adjacency(n_rows);L_columns = chain_adjacency(n_columns);
A = [0,1;1,0];B = [1,0;0,0];C = [0,0;0,1];

adj = kron(kron(kron(I_rows, I_columns), A), O_nodes) + ...
    kron(kron(kron(L_rows, I_columns), B), I_nodes) + ...
    kron(kron(kron(I_rows, L_columns), C), I_nodes);
